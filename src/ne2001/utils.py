"Some utility methods"
from __future__ import division

import numpy as np
from astropy.coordinates import Angle
from astropy.units import Quantity
from numpy import cos
from numpy import pi
from numpy import sin


def rad3d2(xyz):
    return xyz[0]**2 + xyz[1]**2 + xyz[-1]**2


def rad2d2(xyz):
    return xyz[0]**2 + xyz[1]**2


def parse_DM(in_DM):
    """ Convert, as necessary, DM into float
    """
    if isinstance(in_DM, float):
        DM = in_DM
    elif isinstance(in_DM, Quantity):
        DM = in_DM.to('pc/cm**3').value
    else:
        raise IOError("Bad format for input DM")
    # Return
    return DM


def parse_lbd(in_l, in_b, in_d):
    """ Convert, as necessary, l,b,d into floats
    Parameters
    ----------
    in_l : float or Angle
      Galactic longitude; assumed deg if unitless
    in_b : float
      Galactic latitude; assumed deg if unitless
    in_d : float or Quantity
      Distance to source; assumed kpc if unitless

    Returns
    -------
    l : float
    b : float
    d : float
      Distance in kpc

    """
    # l
    if isinstance(in_l, float):
        l = in_l
    elif isinstance(in_l, int):
        l = float(in_l)
    elif isinstance(in_l, (Angle, Quantity)):
        l = in_l.value
    else:
        raise IOError("Bad format for input Galactic longitude")
    # b
    if isinstance(in_b, float):
        b = in_b
    elif isinstance(in_b, int):
        b = float(in_b)
    elif isinstance(in_b, (Angle, Quantity)):
        b = in_b.value
    else:
        raise IOError("Bad format for input Galactic latitude")
    # d
    if isinstance(in_d, float):
        d = in_d
    elif isinstance(in_d, int):
        d = float(in_d)
    elif isinstance(in_d, (Quantity)):
        d = in_d.to('kpc').value
    else:
        raise IOError("Bad format for input distance")
    # Return
    return l, b, d


def galactic_to_galactocentric(l, b, distance, xyz_sun):
    slc = sin(l/180*pi)
    clc = cos(l/180*pi)
    sbc = sin(b/180*pi)
    cbc = cos(b/180*pi)
    rgalc = distance*cbc
    xc = xyz_sun[0] + rgalc*slc
    yc = xyz_sun[1] - rgalc*clc
    zc = xyz_sun[-1] + distance*sbc
    return np.array([xc, yc, zc])


def lzproperty(attribute):
    """
    Lazy property: evaluate property only once
    """
    save_att = '_' + attribute.__name__

    @property
    def _get(self):
        try:
            return getattr(self, save_att)
        except AttributeError:
            setattr(self, save_att, attribute(self))
        return getattr(self, save_att)
    return _get


def rotation(theta, axis=-1):
    """
    Return a rotation matrix around axis
    0:x, 1:y, 2:z
    """
    ct = cos(theta)
    st = sin(theta)

    if axis in (0, -3):
        return np.array([[1, 0, 0],
                         [0, ct, st],
                         [0, -st, ct]])

    if axis in (1, -2):
        return np.array([[ct, 0, st],
                         [0, 1, 0],
                         [-st, 0, ct]])

    if axis in (2, -1):
        return np.array([[ct, st, 0],
                         [-st, ct, 0],
                         [0, 0, 1]])


def rad3d2(xyz):
    return xyz[0]**2 + xyz[1]**2 + xyz[-1]**2


def rad2d2(xyz):
    return xyz[0]**2 + xyz[1]**2


def matmul(a, b):
    try:
        return a.__matmul__(b)
    except AttributeError:
        return np.matmul(a, b)
