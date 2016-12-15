
import cProfile
import os
import pstats
import sys
import time

import numpy as np
from numpy.random import rand
from numpy.random import randint

from ne2001 import density

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../../src/')

PARAMS = density.PARAMS
density.set_xyz_sun(np.array([0, 8.5, 0]))


if __name__ == '__main__':
    tol = 1e-3
    ne = density.ElectronDensity(**PARAMS)
    l, b, d = -2, 12, 1
    DM = 23.98557
    start = time.time()
    DMc = ne.DM(l, b, d)
    t1 = time.time() - start
    start = time.time()
    DMc = ne.DM(l, b, d)
    t2 = time.time() - start
    assert abs(DMc - DM)/DM < tol
    cProfile.run('ne.DM(l, b, d)', 'restats')
    print(t1,t2)
    p = pstats.Stats('restats')
    p.strip_dirs().sort_stats('cumulative').print_stats(50)
