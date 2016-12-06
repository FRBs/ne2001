
import os

import numpy as np
from click.testing import CliRunner
from numpy.random import rand
from numpy.random import randint

from ne2001 import density
from ne2001.cli import main


def test_main():
    runner = CliRunner()
    result = runner.invoke(main, [])

    assert result.output == '()\n'
    assert result.exit_code == 0


def test_density():
    xyz = (1-2*rand(3, 100)) * 20

    ne_disk1 = density.ne_thick_disk(xyz, 0.03, 17.5, 0.97, 8.5)
    assert len(ne_disk1) == 100
    assert all(ne_disk1 >= 0)

    ne_disk2 = density.ne_thin_disk(xyz, 0.08, 3.8, 0.15)
    assert len(ne_disk2) == 100
    assert all(ne_disk2 >= 0)

    ne_gc = density.ne_gc(xyz, 10., 0.145, 0.026,
                          np.array([-0.01, 0.0, -0.029]))
    assert len(ne_gc) == 100
    assert all(ne_gc >= 0)


def test_geometry():
    abc = rand(3)*20
    xyz = (1-2*rand(3, 100)) * 20
    theta = rand()*2*np.pi
    xyz_center = rand(3)*20
    assert density.in_ellisoid(xyz_center, xyz_center, abc, theta)
    assert len(density.in_ellisoid(xyz, xyz_center, abc, theta)) == 100

    abc = rand(3)*20
    xyz = (1-2*rand(3, 100)) * 20
    theta = rand()*2*np.pi
    xyz_center = rand(3)*20
    assert len(density.in_cylinder(xyz, xyz_center, abc, theta)) == 100

    assert all(density.in_half_sphere(xyz, xyz_center, 10)[xyz[-1] < 0] == 0)


def test_clumps():
    clumps_file = os.path.join(os.path.split(density.__file__)[0], "data",
                               "neclumpN.NE2001.dat")

    clumps = density.Clumps(clumps_file)

    clumps = density.Clumps()
    xyz = (clumps.xyz.T[randint(0, clumps.gl.size, 100)].T +
           (1-2*rand(3, 100))*0.01)


    ne_clumps = clumps.ne_clumps(xyz)
    assert len(ne_clumps) == 100
    ix = ne_clumps.argmax()
    assert ne_clumps[ix] == clumps.ne_clumps(xyz[:,ix])


def test_void():
    voids_file = os.path.join(os.path.split(density.__file__)[0], "data",
                              "nevoidN.NE2001.dat")


    voids = density.Voids(voids_file)

    voids = density.Voids()

    xyz = (voids.xyz.T[randint(0, voids.gl.size, 100)].T +
           (1-2*rand(3, 100))*0.01)

    ne_voids = voids.ne_voids(xyz)
    assert len(ne_voids) == 100

    ix = ne_voids.argmax()
    assert ne_voids[ix] == voids.ne_voids(xyz[:,ix])


def test_local_ism():
    xyz = (1-2*rand(3, 100)) * 20
    ne, w = density.ne_local_ism(xyz, density.ldr_params,
                                 density.lsb_params,
                                 density.lhb_params,
                                 density.loop_params)
    assert all(ne >= 0)
    assert len(ne) == 100
    assert len(ne) == len(w)
