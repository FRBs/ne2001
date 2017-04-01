""" Tests on utils module
"""

import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import Angle

from ne2001 import utils


def test_parse_lbd():
    """ Parse lbd """
    # Simple floats
    in_l,in_b,in_d = 1., 1., 50.
    l,b,d = utils.parse_lbd(in_l, in_b, in_d)
    # Test
    for xx in [l,b,d]:
        assert isinstance(xx,float)
    # Angles
    in_l,in_b,in_d = Angle(1.*u.deg), Angle(1.*u.deg), 50.
    l,b,d = utils.parse_lbd(in_l, in_b, in_d)
    assert np.isclose(in_l.value, l)
    # Distance
    in_l,in_b,in_d = Angle(1.*u.deg), Angle(1.*u.deg), 50.*u.kpc
    l,b,d = utils.parse_lbd(in_l, in_b, in_d)
    assert np.isclose(in_d.value, d)
    # Again
    in_l,in_b,in_d = Angle(1.*u.deg), Angle(1.*u.deg), 500.*u.pc
    l,b,d = utils.parse_lbd(in_l, in_b, in_d)
    assert np.isclose(d, 0.5)
    # Arrays
    in_l,in_b,in_d = [1]*3*u.deg, [1]*3*u.deg, [500]*3*u.pc
    l,b,d = utils.parse_lbd(in_l, in_b, in_d)
    assert l.size == b.size == d.size == len(in_l)

    # Bad input
    in_l,in_b,in_d = Angle(1.*u.deg), Angle(1.*u.deg), 500.*u.s
    with pytest.raises(IOError):
        l,b,d = utils.parse_lbd(in_l, in_b, in_d)
    in_l,in_b,in_d = Angle(1.*u.deg), Angle(1.*u.deg), 'abc'
    with pytest.raises(IOError):
        l,b,d = utils.parse_lbd(in_l, in_b, in_d)



def test_parse_DM():
    """ Parse lbd """
    # Simple floats
    in_DM = 20.
    DM = utils.parse_DM(in_DM)
    assert isinstance(DM, float)
    # Quantity
    in_DM = 20. * u.pc / u.cm**3
    DM = utils.parse_DM(in_DM)
    assert np.isclose(DM, in_DM.value)
    # Array
    in_DM = [1]*10*u.pc / u.cm**3
    DM = utils.parse_DM(in_DM)
    assert DM.size == len(DM)
    # Bad inputs
    with pytest.raises(IOError):
        DM = utils.parse_DM('abc')
    with pytest.raises(IOError):
        DM = utils.parse_DM(1*u.s)
