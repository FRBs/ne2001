""" Test I/O algorithms """

import os

import numpy as np
from scipy import integrate

from ne2001 import density
from ne2001 import ne_io


def test_params():
    params = ne_io.Params()
    assert isinstance(params, dict)
    # Thick disk
    assert 'thick_disk' in params
    assert 'thin_disk' in params

def test_galparam():
    gal_param = ne_io.read_galparam()
    assert isinstance(gal_param, dict)
    assert 'thick_disk' in gal_param
    assert 'thin_disk' in gal_param
