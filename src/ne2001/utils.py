"Some utility methods"
from __future__ import division

import numpy as np
from numpy import cos
from numpy import pi
from numpy import sin


def parse_units(val, units, name, dtype=float):
    """
    Convert, as necessary, `val` to `units` as `dtype`
    """

    try:
        return val.to(units).value
    except AttributeError:
        pass
    except Exception as inst:
        raise IOError("Bad format for input {}. ({})".
                      format(name, inst))
    try:
        parsed_val = np.array(val).astype(dtype)
        if parsed_val.size == 1:
            return parsed_val.flatten()[0]
        return parsed_val
    except Exception as inst:
        raise IOError("Bad format for input {}. ({})".
                      format(name, inst))


def parse_DM(DM):
    """
    Convert, as necessary, DM into float
    """
    return parse_units(DM, 'pc/cm**3', 'DM')


def parse_lbd(gal_l, gal_b, distance):
    """
    Convert, as necessary, l,b,d into floats

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
    gal_l : float
      Galactic longitude in deg
    gal_b : float
      Galactic latitude in deg
    distance : float
      Distance in kpc

    """
    l = parse_units(gal_l, 'deg', 'Galactic longitude')
    b = parse_units(gal_b, 'deg', 'Galactic latitude')
    d = parse_units(distance, 'kpc', 'distance')
    return l, b, d


def galactic_to_galactocentric(l, b, distance, xyz_sun):
    """ Convert galactic coordiantes to galactocentric
    Parameters
    ----------
    l : float
      latitude
    b : float
      longitude
    distance : float
      kpc
    xyz_sun : ndarray
      positions of the Sun in kpc

    Returns
    -------
    xyz_c : ndarray
      x,y,z positions along the sightline

    """
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
