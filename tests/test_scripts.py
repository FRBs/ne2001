""" Tests on scripts
"""

import numpy as np
import pytest

from ne2001.scripts import simple_lb


def test_simple_lbd():
    # Simple l,b
    pargs = simple_lb.parser(['1','1'])
    simple_lb.main(pargs)
    # Add d
    pargs = simple_lb.parser(['1','1','-d 50'])
    simple_lb.main(pargs)
