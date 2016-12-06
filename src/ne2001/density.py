"Free electron density model"
import os

import numpy as np
from astropy.table import Table
from numpy import cos
from numpy import cosh
from numpy import exp
from numpy import pi
from numpy import sin
from numpy import sqrt
from numpy import tan


# import astropy.units as us
# from astropy.coordinates import SkyCoord

# Configuration
# TODO: use to config file
# input parameters for large-scale components of NE2001 30 June '02
# flags = {'wg1': 1,
#          'wg2': 1,
#          'wga': 1,
#          'wggc': 1,
#          'wglism': 1,
#          'wgcN': 1,
#          'wgvN': 1}

# solar_params = {'Rsun': 8.3}

# thick_disk_params = {'n1h1': 0.033,
#                      'h1': 0.97,
#                      'A1': 17.5,
#                      'F1': 0.18}

# thin_disk_params = {'n2': 0.08,
#                     'h2': 0.15,
#                     'A2': 3.8,
#                     'F2': 120}

# galactic_center_params = {'xgc': -0.01,
#                           'ygc': 0.0,
#                           'zgc': -0.020,
#                           'rgc': 0.145,
#                           'hgc': 0.026,
#                           'negc0': 10.0,
#                           'Fgc0': 0.6e5}


# spiral_arms_params = {'na': 0.028,
#                       'ha': 0.23,
#                       'wa': 0.65,
#                       'Aa': 10.5,
#                       'Fa': 5,
#                       'narm1': 0.5,
#                       'narm2': 1.2,
#                       'narm3': 1.3,
#                       'narm4': 1.0,
#                       'narm5': 0.25,
#                       'warm1': 1.0,
#                       'warm2': 1.5,
#                       'warm3': 1.0,
#                       'warm4': 0.8,
#                       'warm5': 1.0,
#                       'harm1': 1.0,
#                       'harm2': 0.8,
#                       'harm3': 1.3,
#                       'harm4': 1.5,
#                       'harm5': 1.0,
#                       'farm1': 1.1,
#                       'farm2': 0.3,
#                       'farm3': 0.4,
#                       'farm4': 1.5,
#                       'farm5': 0.3}


ldr_params = {'abc': np.array([1.50, .750, .50]),
              'center': np.array([1.36, 8.06, 0.0]),
              'theta': -24.2*pi/180,
              'ne': 0.012,
              'F': 0.1}

lsb_params = {'abc': np.array([1.050, .4250, .3250]),
              'center': np.array([-0.75, 9.0, -0.05]),
              'theta': 139.*pi/180,
              'ne': 0.016,
              'F':  0.01}

lhb_params = {'abc': np.array([.0850, .1000, .330]),
              'center': np.array([0.01, 8.45, 0.17]),
              'theta': 15*pi/180,
              'ne': 0.005,
              'F': 0.01}

loop_params = {'center': np.array([-0.045, 8.40, 0.07]),
               'r': 0.120,
               'dr': 0.060,
               'ne1': 0.0125,
               'ne2': 0.0125,
               'F1': 0.2,
               'F2': 0.01}


def ne_thick_disk(xyz, ne_disk, rdisk, hdisk, r_sun):
    """
    Calculate the contribution of the thick disk to the free electron density
     at x, y, z = `xyz`
    """
    r2d = sqrt(xyz[0]**2 + xyz[1]**2)
    k = pi/2/rdisk
    return ne_disk * (cos(r2d*k)/cos(r_sun*k) /
                      cosh(xyz[-1]/hdisk)**2 *
                      (r2d < rdisk))


def ne_thin_disk(xyz, ne_disk, rdisk, hdisk):
    """
    Calculate the contribution of the thin disk to the free electron density
     at x, y, z = `xyz`
    """
    r2d = sqrt(xyz[0]**2 + xyz[1]**2)
    return ne_disk * (exp(-(r2d - rdisk)**2/1.8**2) /
                      cosh(xyz[-1]/hdisk)**2)  # Why 1.8?


def ne_gc(xyz, ne_gc0, rgc, hgc, xyz_gc):
    """
    Calculate the contribution of the Galactic center to the free
    electron density at x, y, z = `xyz`
    """
    # Here I'm using the expression in the NE2001 code which is inconsistent
    # with Cordes and Lazio 2011 (0207156v3) (See Table 2)
    try:
        xyz = xyz - xyz_gc
    except ValueError:
        xyz = xyz - xyz_gc[:, None]

    r2d = sqrt(xyz[0]**2 + xyz[1]**2)

    # ????
    # Cordes and Lazio 2011 (0207156v3) (Table 2)
    # return ne_gc0*exp(-(r2d/rgc)**2 - (xyz[-1]/hgc)**2)
    # ????

    # Constant ne (form NE2001 code)
    return ne_gc0*((r2d/rgc)**2 + (xyz[-1]/hgc)**2 < 1)*(r2d < rgc)


def ne_local_ism(xyz, ldr_params, lsb_params, lhb_params, loop_params):
    """
    Calculate the contribution of the local ISM to the free
    electron density at x, y, z = `xyz`
    """
    # low density region in Q1
    neldr = ldr_params['ne']*in_ellisoid(xyz, ldr_params['center'],
                                         ldr_params['abc'],
                                         ldr_params['theta'])
    # Local Super Bubble
    nelsb = lsb_params['ne']*in_ellisoid(xyz, lsb_params['center'],
                                         lsb_params['abc'],
                                         lsb_params['theta'])

    # Local Hot Bubble
    nelhb = lhb_params['ne']*in_cylinder(xyz, lhb_params['center'],
                                         lhb_params['abc'],
                                         lhb_params['theta'])
    # Loop I
    irr1 = in_half_sphere(xyz, loop_params['center'], loop_params['r'])
    irr2 = in_half_sphere(xyz, loop_params['center'],
                          loop_params['r'] + loop_params['dr'])
    neloop = loop_params['ne1'] * irr1 + loop_params['ne2'] * irr2*(~irr1)
    wlhb, wloop, wlsb, wldr = (nelhb > 0,
                               neloop > 0,
                               nelsb > 0,
                               neldr > 0)
    ne_lism = ((1 - wlhb) *
               ((1 - wloop) * (wlsb*nelsb + (1-wlsb) * neldr) +
                wloop*neloop) +
               wlhb*nelhb)

    wlism = np.maximum(wloop, np.maximum(wldr, np.maximum(wlsb, wlhb)))
    return ne_lism, wlism


def in_ellisoid(xyz, xyz_center, abc_ellipsoid, theta):
    """
    Test if xyz in the ellipsoid
    Theta in radians
    """
    try:
        xyz = xyz - xyz_center
    except ValueError:
        xyz = xyz - xyz_center[:, None]
        abc_ellipsoid = abc_ellipsoid[:, None]

    rot = rotation(theta, -1)
    xyz = rot.dot(xyz)

    xyz_p = xyz/abc_ellipsoid

    return np.sum(xyz_p**2, axis=0) <= 1


def in_cylinder(xyz, xyz_center, abc_cylinder, theta):
    """
    Test if xyz in the cylinder
    Theta in radians
    """
    try:
        xyz = xyz - xyz_center
    except ValueError:
        xyz = xyz - xyz_center[:, None]
        abc_cylinder = np.vstack([abc_cylinder]*xyz.shape[-1]).T
    xyz[2] -= tan(theta)*xyz[-1]

    abc_cylinder_p = abc_cylinder.copy()
    z_c = (xyz_center[-1] - abc_cylinder[-1])
    izz = (xyz[-1] <= 0)*(xyz[-1] <= z_c)
    abc_cylinder_p[0] = (0.001 +
                         (abc_cylinder[0] - 0.001) *
                         (1 - xyz[-1]/z_c))*izz + abc_cylinder[0]*(~izz)

    xyz_p = xyz/abc_cylinder_p

    return (xyz_p[0]**2 + xyz_p[1]**2 <= 1) * (xyz_p[-1]**2 <= 1)


def in_half_sphere(xyz, xyz_center, r_sphere):
    "Test if `xyz` in the sphere with radius r_sphere  centerd at `xyz_center`"
    try:
        xyz = xyz - xyz_center
    except ValueError:
        xyz = xyz - xyz_center[:, None]
    distance = sqrt(np.sum(xyz**2, axis=0))
    return (distance <= r_sphere)*(xyz[-1] >= 0)


class Clumps(object):
    """
    """

    def __init__(self, clumps_file=None):
        """
        """
        if not clumps_file:
            this_dir, _ = os.path.split(__file__)
            clumps_file = os.path.join(this_dir, "data", "neclumpN.NE2001.dat")
        self._data = Table.read(clumps_file, format='ascii')

    @property
    def use_clump(self):
        """
        """
        return self._data['flag'] == 0

    @property
    def xyz(self):
        """
        """
        try:
            return self._xyz
        except AttributeError:
            self._xyz = self.get_xyz()
        return self._xyz

    @property
    def gl(self):
        """
        Galactic longitude (deg)
        """
        return self._data['l']

    @property
    def gb(self):
        """
        Galactic latitude (deg)
        """
        return self._data['b']

    @property
    def distance(self):
        """
        Distance from the sun (kpc)
        """
        return self._data['dc']

    @property
    def radius(self):
        """
        Radius of the clump (kpc)
        """
        return self._data['rc']

    @property
    def ne(self):
        """
        Electron density of each clump (cm^{-3})
        """
        return self._data['nec']

    @property
    def edge(self):
        """
        Clump edge
        0 => use exponential rolloff out to 5 clump radii
        1 => uniform and truncated at 1/e clump radius
        """
        return self._data['edge']

    def get_xyz(self, rsun=8.5):
        """
        """
        # xyz = SkyCoord(frame="galactic", l=self.gl, b=self.gb,
        #                distance=self.distance,
        #                z_sun = z_sun*us.kpc,
        #                unit="deg, deg, kpc").galactocentric.cartesian.xyz.value
        # return xyz

        slc = sin(self.gl/180*pi)
        clc = cos(self.gl/180*pi)
        sbc = sin(self.gb/180*pi)
        cbc = cos(self.gb/180*pi)
        rgalc = self.distance*cbc
        xc = rgalc*slc
        yc = rsun-rgalc*clc
        zc = self.distance*sbc
        return np.array([xc, yc, zc])

    def clump_factor(self, xyz):
        """
        Clump edge
        0 => use exponential rolloff out to 5 clump radii
        1 => uniform and truncated at 1/e clump radius
        """
        if xyz.ndim == 1:
            xyz = xyz[:, None] - self.xyz
        else:
            xyz = xyz[:, :, None] - self.xyz[:, None, :]

        q2 = (np.sum(xyz**2, axis=0) /
              self.radius**2)
        # NOTE: In the original NE2001 code q2 <= 5 is used instead of q <= 5.
        # TODO: check this
        return (q2 <= 1)*(self.edge == 1) + (q2 <= 5)*(self.edge == 0)*exp(-q2)

    def ne_clumps(self, xyz):
        """
        The contribution of the clumps to the free
        electron density at x, y, z = `xyz`
        """
        return np.sum(self.clump_factor(xyz)*self.ne*self.use_clump, axis=-1)


class Voids(object):
    """
    """

    def __init__(self, voids_file=None):
        """
        """
        if not voids_file:
            this_dir, _ = os.path.split(__file__)
            voids_file = os.path.join(this_dir, "data", "nevoidN.NE2001.dat")
        self._data = Table.read(voids_file, format='ascii')

    @property
    def use_void(self):
        """
        """
        return self._data['flag'] == 0

    @property
    def xyz(self):
        """
        """
        try:
            return self._xyz
        except AttributeError:
            self._xyz = self.get_xyz()
        return self._xyz

    @property
    def gl(self):
        """
        Galactic longitude (deg)
        """
        return self._data['l']

    @property
    def gb(self):
        """
        Galactic latitude (deg)
        """
        return self._data['b']

    @property
    def distance(self):
        """
        Distance from the sun (kpc)
        """
        return self._data['dv']

    @property
    def ellipsoid_abc(self):
        """
        Void axis
        """
        return np.array([self._data['aav'],
                         self._data['bbv'],
                         self._data['ccv']])

    @property
    def rotation_y(self):
        """
        Rotation around the y axis
        """
        return [rotation(theta*pi/180, 1) for theta in self._data['thvy']]

    @property
    def rotation_z(self):
        """
        Rotation around the z axis
        """
        return [rotation(theta*pi/180, -1) for theta in self._data['thvz']]

    @property
    def ne(self):
        """
        Electron density of each void (cm^{-3})
        """
        return self._data['nev']

    @property
    def edge(self):
        """
        Void edge
        0 => use exponential rolloff out to 5 clump radii
        1 => uniform and truncated at 1/e clump radius
        """
        return self._data['edge']

    def get_xyz(self, rsun=8.5):
        """
        """
        # xyz = SkyCoord(frame="galactic", l=self.gl, b=self.gb,
        #                distance=self.distance,
        #                z_sun = z_sun*us.kpc,
        #                unit="deg, deg, kpc").galactocentric.cartesian.xyz.value
        # return xyz

        slc = sin(self.gl/180*pi)
        clc = cos(self.gl/180*pi)
        sbc = sin(self.gb/180*pi)
        cbc = cos(self.gb/180*pi)
        rgalc = self.distance*cbc
        xc = rgalc*slc
        yc = rsun-rgalc*clc
        zc = self.distance*sbc
        return np.array([xc, yc, zc])

    def void_factor(self, xyz):
        """
        Clump edge
        0 => use exponential rolloff out to 5 clump radii
        1 => uniform and truncated at 1/e clump radius
        """
        if xyz.ndim == 1:
            xyz = xyz[:, None] - self.xyz
            ellipsoid_abc = self.ellipsoid_abc
        else:
            xyz = xyz[:, :, None] - self.xyz[:, None, :]
            ellipsoid_abc = self.ellipsoid_abc[:, None, :]

        xyz = np.array([Rz.dot(Ry).dot(XYZ.T).T
                        for Rz, Ry, XYZ in
                        zip(self.rotation_z, self.rotation_y, xyz.T)]).T

        q2 = np.sum(xyz**2 / ellipsoid_abc**2, axis=0)
        # NOTE: In the original NE2001 code q2 <= 5 is used instead of q <= 5.
        # TODO: check this
        return (q2 <= 1)*(self.edge == 1) + (q2 <= 5)*(self.edge == 0)*exp(-q2)

    def ne_voids(self, xyz):
        """
        The contribution of the clumps to the free
        electron density at x, y, z = `xyz`
        """
        return np.sum(self.void_factor(xyz)*self.ne*self.use_void, axis=-1)


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
