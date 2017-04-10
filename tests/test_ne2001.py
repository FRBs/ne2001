
import os

import numpy as np
import pytest
from click.testing import CliRunner
from numpy.random import rand
from numpy.random import randint
from numpy.random import seed
from scipy import integrate

from ne2001 import density
from ne2001 import ne_io
from ne2001 import utils
from ne2001.cli import main

PARAMS = ne_io.Params()
density.set_xyz_sun(np.array([0, 8.5, 0]))


def test_main():
    runner = CliRunner()
    result = runner.invoke(main, [])

    assert result.output == '()\n'
    assert result.exit_code == 0


def test_density():
    xyz = (1-2*rand(3, 100)) * 20

    ne_disk1 = density.NEobject(density.thick_disk,
                                **PARAMS['thick_disk']).ne(xyz)
    assert len(ne_disk1) == 100
    assert all(ne_disk1 >= 0)

    ne_disk2 = density.NEobject(density.thin_disk, **PARAMS['thin_disk']).ne(xyz)
    assert ne_disk2.size == 100
    assert all(ne_disk2 >= 0)

    ne_gc = density.NEobject(density.gc, **PARAMS['galactic_center']).ne(xyz)
    assert len(ne_gc) == 100
    assert all(ne_gc >= 0)


def test_geometry():
    ellipsoid = rand(3)*20
    xyz = (1-2*rand(3, 100)) * 20
    theta = rand()*2*np.pi
    center = rand(3)*20
    assert density.in_ellipsoid(center, ellipsoid, theta)(center)
    assert len(density.in_ellipsoid(center, ellipsoid, theta)(xyz)) == 100

    cylinder = rand(3)*20
    xyz = (1-2*rand(3, 100)) * 20
    theta = rand()*2*np.pi
    center = rand(3)*20
    assert len(density.in_cylinder(xyz, center, cylinder, theta)) == 100

    assert all(density.in_half_sphere(xyz, center, 10)[xyz[-1] < 0] == 0)


def test_clumps():
    clumps_file = os.path.join(os.path.split(density.__file__)[0], "data",
                               "neclumpN.NE2001.dat")

    clumps = density.Clumps(clumps_file)

    clumps = density.Clumps()
    xyz = (clumps.xyz.T[randint(0, clumps.gl.size, 100)].T +
           (1-2*rand(3, 100))*0.01)

    ne_clumps = clumps.ne(xyz)
    assert len(ne_clumps) == 100
    ix = ne_clumps.argmax()
    assert ne_clumps[ix] == clumps.ne(xyz[:, ix])


def test_void():
    tol = 1e-3
    voids_file = os.path.join(os.path.split(density.__file__)[0], "data",
                              "nevoidN.NE2001.dat")

    voids = density.Voids(voids_file)

    voids = density.Voids()

    xyz = (voids.xyz.T[randint(0, voids.gl.size, 100)].T +
           (1-2*rand(3, 100))*0.01)

    ne_voids = voids.ne(xyz)
    assert len(ne_voids) == 100

    ix = ne_voids.argmax()
    assert ne_voids[ix] == voids.ne(xyz[:, ix])
    l, b, d = -2, 12, 1
    DM = 0.9308054
    assert abs(voids.DM(l, b, d).value - DM)/DM < tol


def test_local_ism():
    tol = 1e-3
    xyz = (1-2*rand(3, 100)) * 20
    local_ism = density.LocalISM(**PARAMS)

    assert all(local_ism.ne(xyz) >= 0)
    assert len(local_ism.ne(xyz)) == 100
    l, b, d = -2, 12, 1
    DM = 2.453550
    assert (abs(local_ism.DM(l, b, d).value -
                DM)/DM < tol)

def test_DM():
    tol = 1e-3
    d1 = density.NEobject(density.thick_disk, **PARAMS['thick_disk'])
    xyz = np.array([0.1503843, 7.647129, 0.5000018])
    l, b, d = 10, 30, 1
    DM = 32.36372
    assert abs(d1.DM(l, b, d).value - DM) / DM  < tol


    d2 = density.NEobject(density.thin_disk,
                          **PARAMS['thin_disk'])


    l, b, d = 10, 30, 1
    DM = 0.046937
    assert abs(d2.DM(l, b, d).value - DM) / DM  < tol


def test_electron_density_quad():
    tol = 1e-3
    ne = density.ElectronDensity()
    l, b, d = -2, 12, 1
    DM = 23.98557
    assert abs(ne.DM(l, b, d).value - DM)/DM < tol


def test_electron_density_trapz():
    tol = 1e-3
    ne = density.ElectronDensity()
    l, b, d = -2, 12, 1
    DM = 23.98557
    assert abs(ne.DM(l, b, d, integrator=integrate.trapz).value - DM)/DM < tol


def test_dist():
    seed(123)
    for i in range(1):
        tol = 0.1
        ne = density.ElectronDensity()
        l = rand()*360
        b = np.arccos(1 - 2*rand())/np.pi*180
        d = rand()*5
        DM = ne.DM(l, b, d)
        d_DM = ne.dist(l, b, DM)
        err = abs(d_DM.value - d)/d
        print(err, l, b, d, d_DM)
        assert err < tol, (l, b, d)
