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


class ClassOperation(object):
    """
    Define arithmetic operation on Class
    """
    def __add__(self, other):
        return ClassOperationResult('__add__', self, other)

    def __sub__(self, other):
        return ClassOperationResult('__sub__', self, other)

    def __mul__(self, other):
        return ClassOperationResult('__mul__', self, other)

    def __gt__(self, other):
        return ClassOperationResult('__gt__', self, other)


class ClassOperationResult(ClassOperation):
    """
    Result of Class Operation
    """
    def __init__(self, operation, cls1, cls2):
        """
        """
        self.cls1 = cls1
        self.cls2 = cls2
        self._operation = operation

    def __getattr__(self, attr):
        attr1 = getattr(self.cls1, attr)
        attr2 = getattr(self.cls2, attr)

        try:
            return getattr(attr1, self._operation)(attr2)

        except AttributeError:
            return (lambda *args, **kwargs:
                    getattr(attr1(*args), self._operation)(attr2(*args)))


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


def galactic_to_galactocentric(l, b, distance, rsun):
    slc = sin(l/180*pi)
    clc = cos(l/180*pi)
    sbc = sin(b/180*pi)
    cbc = cos(b/180*pi)
    rgalc = distance*cbc
    xc = rgalc*slc
    yc = rsun-rgalc*clc
    zc = distance*sbc
    return np.array([xc, yc, zc])
