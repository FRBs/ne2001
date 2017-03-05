""" Module for I/O
"""
from __future__ import absolute_import, division, print_function

import numpy as np
import json
import pdb

import ne2001
data_path = ne2001.__path__[0]+'/data/'


def read_params(ifile='ne2001_params.json'):
    """ Read parameter file"""
    # Read
    with open(data_path+ifile, 'rt') as fh:
        PARAMS = json.load(fh)
    # Add
    PARAMS['spiral_arms']['adict']  = init_spiral_arms()
    PARAMS['spiral_arms']['gal_param'] = read_galparam()
    # Recast?
    return PARAMS


def read_galparam(ifile='gal_param.json'):
    """ Read Galaxy parameters
    Parameters
    ----------
    ifile : str, optional

    Returns
    -------
    gal_param : dict

    """
    # Read
    with open(data_path+ifile, 'rt') as fh:
        galparam_dict = json.load(fh)
    # Thick disk
    thick_dict = dict(e_density=galparam_dict['n1h1']/galparam_dict['h1'],
                      height=galparam_dict['h1'],
                      radius=galparam_dict['A1'],
                      F=galparam_dict['F1'],
                      )
    # Thin disk
    thin_dict = dict(e_density=galparam_dict['n2'],
                      height=galparam_dict['h2'],
                      radius=galparam_dict['A2'],
                      F=galparam_dict['F2'],
                      )
    # Return
    return galparam_dict


def read_gc(ifile='ne_gc.json'):
    """ Read Galactic Center parameters
    Returns
    -------
    gc_dict : dict
      dict of parameters

    """
    # Read
    with open(data_path+ifile, 'rt') as fh:
        gc_dict = json.load(fh)
    # Return
    return gc_dict


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
    with open(data_path+ifile, 'rt') as fh:
        lism_dict = json.load(fh)
    # Return
    return lism_dict


def init_spiral_arms():

    from astropy.table import Table
    from scipy.interpolate import CubicSpline
    armsinp= data_path + 'ne_arms_log_mod.inp'
    logarms= data_path + 'log_arms.out'

    narms=5
    #integer armmap(5)		! for remapping from Wainscoat
    #data armmap/1, 3, 4, 2, 5/	! order to TC93 order, which is
    #                ! from GC outwards toward Sun.
    armmap = [1,3,4,2,5]
    NNj = [20, 20, 20, 20, 20]
    narmpoints=500
    ncoord=2
    NNmax=20
    rad = 180/np.pi
    # Arms
    arms_tbl = Table.read(armsinp, format='ascii') # a, rmin, thmin, extent
    assert len(arms_tbl) == narms

    r1 = np.zeros((NNmax, narms))
    th1 = np.zeros((NNmax, narms))
    kmax = np.zeros(narms).astype(int)
    arm = np.zeros((narms, narmpoints, ncoord))

    for j, row in enumerate(arms_tbl):
        th1[0:NNj[j],j] = row['thmin'] + np.arange(NNj[j])*row['extent']/(NNj[j]-1.) 	#! rad
        r1[:,j] = row['rmin'] * np.exp((th1[:,j]-row['thmin'])/row['a'])
        th1[:,j] *= rad # ! deg
        #c *** begin sculpting spiral arm 2 == TC arm 3***
        if armmap[j] == 3:
            cut1 = (th1[:,j] > 370.) & (th1[:,j] <= 410.)
            r1[cut1,j] *= (1. + 0.04* np.cos((th1[cut1,j]-390.)*180./(40.*rad)))
                #c    .           (1. + 0.01*cos((th1(n,j)-390.)*180./(40.*rad)))
            cut2 = (th1[:,j] > 315.) & (th1[:,j] <= 370.)
            r1[cut2,j] *= (1. - 0.07* np.cos((th1[cut2,j]-345.)*180./(55.*rad)))
                #c    .           (1.0 - 0.08*cos((th1(n,j)-345.)*180./(55.*rad)))
            cut3 = (th1[:,j] > 180.) & (th1[:,j] <= 315.)
            r1[cut3,j] *= (1 + 0.16* np.cos((th1[cut3,j]-260.)*180./(135.*rad)))
                    # (1 + 0.13* np.cos((th1[cut3,j]-260.)*180./(135.*rad)))
        #c *** begin sculpting spiral arm 4 == TC arm 2***
        if armmap[j] == 2:
            cut1 = (th1[:,j] > 290.) & (th1[:,j] <= 395.)
            r1[cut1,j] *= (1. - 0.11* np.cos((th1[cut1,j]-350.)*180./(105.*rad)))
        #c *** end arm sculpting ***



    """
    open(11,file=logarms, status='unknown')
        write(11,*) 'arm  n   xa     ya'
    """
    #do 21 j=1,narms
    for j in range(narms):
        dth = 5.0/r1[0,j]  # Python indexing
        th = th1[0,j]-0.999*dth
        # Generate spline
        cspline = CubicSpline(th1[:NNj[j],j],r1[:NNj[j],j])
        #call cspline(th1(1,j),r1(1,j),-NNj(j),th,r)
        #for k in range(narmpoints):
          #do 10 k=1,narmpoints-1
        th = th + dth * np.arange(narmpoints)
        gd_th = np.where(th <= th1[NNj[j]-1, j])[0]
        kmax[j] = np.max(gd_th) + 1  # Python indexing (we will use arange)
        r = cspline(th[gd_th])
        # x,y of each arm
        arm[j,gd_th,0] = -r*np.sin(th[gd_th]/rad)  # Python indexing
        arm[j,gd_th,1] = r*np.cos(th[gd_th]/rad)

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
