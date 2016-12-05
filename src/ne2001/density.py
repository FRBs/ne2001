"Free electron density model"
import numpy as np
from numpy import cos
from numpy import cosh
from numpy import exp
from numpy import pi
from numpy import sin
from numpy import sqrt


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

def ne_thick_disk(xyz, ne_disk, hdisk, rdisk, r_sun):
    """
    Calculate the contribution of the thick disk to the free electron density
     at x, y, z = `xyz`
    """
    r2d = sqrt(xyz[0]**2 + xyz[1]**2)
    k = pi/2/rdisk
    return ne_disk * (cos(r2d*k)/cos(r_sun*k) /
                      cosh(xyz[-1]/hdisk)**2 *
                      (r2d < rdisk))


def ne_thin_disk(xyz, ne_disk, hdisk, rdisk):
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
    xyz = (xyz_gc - xyz)
    r2d = sqrt(xyz[0]**2 + xyz[1]**2)

    # ????
    # Cordes and Lazio 2011 (0207156v3) (Table 2)
    # return ne_gc0*exp(-(r2d/rgc)**2 - (xyz[-1]/hgc)**2)
    # ????

    # Constant ne (form NE2001 code)
    return ne_gc0*((r2d/rgc)**2 + (xyz[-1]/hgc)**2 < 1)*(r2d < rgc)


def in_ellisoid(xyz, xyz_center, abc_ellipsoid, theta):
    """
    Test if xyz in the ellipsoid
    """
    xyz = xyz - xyz_center
    rot = np.array([[cos(theta), sin(theta)],
                    [-sin(theta), cos(theta)]])
    xyz[:2] = xyz @ rot

    xyz_p = xyz/abc_ellipsoid

    return np.sum(xyz_p**2, axis=0) < 1
