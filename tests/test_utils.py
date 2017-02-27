""" Tests on utils module
"""

import numpy as np
import pytest

from ne2001 import utils

from astropy.coordinates import Angle
from astropy import units as u


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


