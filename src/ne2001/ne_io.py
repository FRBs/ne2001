""" Module for I/O
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import json
import os
from builtins import super

import numpy as np
from astropy.table import Table
from scipy.interpolate import CubicSpline

from . import __path__

DATA_PATH = os.path.join(__path__[0], 'data')


def numpify_dict(d):
    """
    Recursively make lists in a dictionary into numpy array
    """
    def numpify(d):
        for k, v in d.items():
            if isinstance(v, list):
                d[k] = np.array(v)
            elif isinstance(v, dict):
                numpify(v)
    new_dict = d.copy()
    numpify(new_dict)
    return new_dict


class Params(dict):
    """
    Input parameters
    """

    def __init__(self, ifile='ne2001_params.json', path=None, **new_params):
        """
        """
        if path is None:
            path = DATA_PATH
        self.path = path
        self.ifile = ifile
        try:
            params = numpify_dict(parse_json(os.path.join(self.path,
                                                          self.ifile)))
            params['spiral_arms']['adict'] = init_spiral_arms()
        except IOError:
            params = {}
        params.update(new_params)
        super().__init__(params)


def parse_json(json_file):
    "Parse json file"
    with open(json_file, 'rt') as json_data:
        data = json.load(json_data)
    return data


def read_galparam(ifile='gal_param.json'):
    """
    Read Galaxy parameters

    Parameters
    ----------
    ifile : str, optional

    Returns
    -------
    gal_param : dict

    """
    old_param = parse_json(os.path.join(DATA_PATH, ifile))
    gal_param = {}

    gal_param['thick_disk'] = dict(e_density=(old_param['n1h1'] /
                                              old_param['h1']),
                                   height=old_param['h1'],
                                   radius=old_param['A1'],
                                   F=old_param['F1'])

    gal_param['thin_disk'] = dict(e_density=old_param['n2'],
                                  height=old_param['h2'],
                                  radius=old_param['A2'],
                                  F=old_param['F2'])

    return gal_param


def read_gc(ifile='ne_gc.json'):
    """ Read Galactic Center parameters
    Returns
    -------
    gc_param : dict
      dict of parameters

    """
    old_param = parse_json(os.path.join(DATA_PATH, ifile))
    gc_param = {}

    gc_param['galactic_center'] = dict(e_density=old_param['negc0'],
                                       center=tuple(old_param['centroid'].
                                                    values()),
                                       F=old_param['Fgc0'],
                                       height=old_param['hgc'],
                                       radius=old_param['rgc'])

    return gc_param


def read_lism(ifile='ne_lism.json'):
    """
    Parameters
    ----------
    ifile : str, optional

    Returns
    -------
    lism_dict : dict

    """
    # Read
    with open(os.path.join(DATA_PATH, ifile), 'rt') as fh:
        lism_dict = json.load(fh)
    # Return
    return lism_dict


def init_spiral_arms(ifile='ne_arms_log_mod.inp'):
    armsinp = os.path.join(DATA_PATH, ifile)
    # logarms = DATA_PATH + 'log_arms.out'

    narms = 5
    # integer armmap(5)		! for remapping from Wainscoat
    # data armmap/1, 3, 4, 2, 5/	! order to TC93 order, which is
    #                ! from GC outwards toward Sun.
    armmap = [1, 3, 4, 2, 5]
    NNj = [20, 20, 20, 20, 20]
    narmpoints = 500
    ncoord = 2
    NNmax = 20
    rad = 180/np.pi
    # Arms
    arms_tbl = Table.read(armsinp, format='ascii')  # a, rmin, thmin, extent
    assert len(arms_tbl) == narms

    r1 = np.zeros((NNmax, narms))
    th1 = np.zeros((NNmax, narms))
    kmax = np.zeros(narms).astype(int)
    arm = np.zeros((narms, narmpoints, ncoord))

    for j, row in enumerate(arms_tbl):
        th1[0:NNj[j], j] = (row['thmin'] +
                            np.arange(NNj[j])*row['extent']/(NNj[j]-1.))  # rad
        r1[:, j] = row['rmin'] * np.exp((th1[:, j]-row['thmin'])/row['a'])
        th1[:, j] *= rad  # ! deg
        # c *** begin sculpting spiral arm 2 == TC arm 3***
        if armmap[j] == 3:
            cut1 = (th1[:, j] > 370.) & (th1[:, j] <= 410.)
            r1[cut1, j] *= (1. + 0.04 * np.cos((th1[cut1, j]-390.)*180 /
                                               (40.*rad)))
            # c . (1. + 0.01*cos((th1(n,j)-390.)*180./(40.*rad)))
            cut2 = (th1[:, j] > 315.) & (th1[:, j] <= 370.)
            r1[cut2, j] *= (1. - 0.07 * np.cos((th1[cut2, j]-345.)*180 /
                                               (55.*rad)))
            # c   .   (1.0 - 0.08*cos((th1(n,j)-345.)*180./(55.*rad)))
            cut3 = (th1[:, j] > 180.) & (th1[:, j] <= 315.)
            r1[cut3, j] *= (1 + 0.16 * np.cos((th1[cut3, j]-260.)*180 /
                                              (135.*rad)))
            # (1 + 0.13* np.cos((th1[cut3,j]-260.)*180./(135.*rad)))
        # c *** begin sculpting spiral arm 4 == TC arm 2***
        if armmap[j] == 2:
            cut1 = (th1[:, j] > 290.) & (th1[:, j] <= 395.)
            r1[cut1, j] *= (1. - 0.11 * np.cos((th1[cut1, j]-350.)*180 /
                                               (105.*rad)))
        # c *** end arm sculpting ***

    """
    open(11,file=logarms, status='unknown')
        write(11,*) 'arm  n   xa     ya'
    """
    # do 21 j=1,narms
    for j in range(narms):
        dth = 5.0/r1[0, j]  # Python indexing
        th = th1[0, j]-0.999*dth
        # Generate spline
        cspline = CubicSpline(th1[:NNj[j], j], r1[:NNj[j], j])
        # call cspline(th1(1,j),r1(1,j),-NNj(j),th,r)
        # for k in range(narmpoints):
        # do 10 k=1,narmpoints-1
        th = th + dth * np.arange(narmpoints)
        gd_th = np.where(th <= th1[NNj[j]-1, j])[0]
        kmax[j] = np.max(gd_th) + 1  # Python indexing (we will use arange)
        r = cspline(th[gd_th])
        # x,y of each arm
        arm[j, gd_th, 0] = -r*np.sin(th[gd_th]/rad)  # Python indexing
        arm[j, gd_th, 1] = r*np.cos(th[gd_th]/rad)

    # Wrap into a dict
    arms_dict = {}
    arms_dict['table'] = arms_tbl
    arms_dict['r1'] = r1
    arms_dict['th1'] = r1
    arms_dict['kmax'] = kmax
    arms_dict['narms'] = narms
    arms_dict['narmpoints'] = narmpoints
    arms_dict['armmap'] = armmap
    arms_dict['arm'] = arm
    return arms_dict
