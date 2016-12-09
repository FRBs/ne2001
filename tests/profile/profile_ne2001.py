
import cProfile
import os
import pstats
import time

import numpy as np
from click.testing import CliRunner
from numpy.random import rand
from numpy.random import randint

from ne2001 import density
from ne2001.cli import main

PARAMS = density.PARAMS
r_sun = 8.5


if __name__ == '__main__':
    tol = 1e-3
    ne = density.ElectronDensity(**PARAMS)
    xyz = np.array([-3.4153607E-02,   7.521969,      0.2080137])
    DM = 23.98557
    start = time.time()
    DMc = ne.DM(xyz)
    t1 = time.time() - start
    start = time.time()
    DMc = ne.DM(xyz)
    t2 = time.time() - start
#    assert abs(DMc - DM)/DM < tol
    cProfile.run('ne.DM(xyz)', 'restats')
    print(t1,t2)
    p = pstats.Stats('restats')
    p.strip_dirs().sort_stats('cumulative').print_stats(50)
