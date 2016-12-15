"Some utility methods"
from __future__ import division
import numpy as np
from numpy import cos
from numpy import pi
from numpy import sin


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


def galactic_to_galactocentric(l, b, distance, xyz_sun=[0, 8.5, 0]):
    slc = sin(l/180*pi)
    clc = cos(l/180*pi)
    sbc = sin(b/180*pi)
    cbc = cos(b/180*pi)
    rgalc = distance*cbc
    xc = xyz_sun[0] + rgalc*slc
    yc = xyz_sun[1] - rgalc*clc
    zc = xyz_sun[-1] + distance*sbc
    return np.array([xc, yc, zc])
