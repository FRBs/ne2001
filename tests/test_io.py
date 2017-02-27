""" Test I/O algorithms """

import os

import numpy as np
from scipy import integrate

from ne2001 import density
from ne2001 import io as ne_io

def test_params():
    PARAMS = ne_io.read_params()
    assert isinstance(PARAMS, dict)
    # Thick disk
    assert 'thick_disk' in PARAMS.keys()

def test_galparam():
    gal_param = ne_io.read_galparam()
    assert isinstance(gal_param, dict)

    for key in ['A1','A2','harm1','h1']:
        assert key in gal_param.keys()

