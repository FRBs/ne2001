""" Density for spiral arms
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
from scipy.interpolate import CubicSpline

from .utils import rad2d2


def ne_spiral_arm(xyz, Aa, wa, ha, farms, harms, narms, warms, adict):
    """
    Parameters
    ----------
    xyz
    gal_param
    adict

    Returns
    -------
    nea : ndarray or float
    whicharm : ndarray or int

    c-----------------------------------------------------------------------
    c  Spiral arms are defined as logarithmic spirals using the
    c    parameterization in Wainscoat et al. 1992, ApJS, 83, 111-146.
    c  But arms are modified selectively at various places to distort them
    c    as needed (08 Aug 2000).
    c  Note that arm numbering follows that of TC93 for the four large arms
    c (after remapping).
    c  The local spiral arm is number 5.
    c  06 Apr 02:   removed TC type modifications of arms 2,3 (fac calculations)
    c  		and replaced with new versions.  Data for these are hard wired.
    """
    x, y, z = xyz[0], xyz[1], xyz[-1]
    if isinstance(x, float):
        x = np.array([x])
        y = np.array([y])
        z = np.array([z])
        flg_float = True
    else:
        flg_float = False

    # see get_parameters for definitions of narm, warm, harm.
    # narmsmax = 5
    # common/armfactors/
    # .     harm(narmsmax),narm(narmsmax),warm(narmsmax),farm(narmsmax)

    #   parameter(rad=57.29577 95130 823)
    rad = 180/np.pi
    # ks = 3
    # NN = 7
    nfine = 1000

    # rr
    rr = rad2d2(xyz)
    if flg_float:
        rr = np.array([rr])

    #
    # Get spiral arm component:  30 do loop finds a coarse minimum distance
    # from line of sight to arm; 40 do loop finds a fine minimum distance
    # from line of sight to arm; line 35 ensures that arm limits are not
    # exceeded; linear interpolation beginning at line 41 finds the
    # minimum distance from line of sight to arm on a finer scale than gridding
    # of arms allows (TJL)

    # Init
    whicharm = np.zeros_like(x).astype(int)
    nea = np.zeros_like(x)

    # thxy
    thxy = np.arctan2(-x, y) * rad  # measured ccw from +y axis (different from tc93 theta)
    neg_th = thxy < 0.
    thxy[neg_th] += 360.
    # Cut on values near the disk
    icutz = np.where(np.abs(z/ha) < 10.)[0]
    if len(icutz) > 0:
        cutx = x[icutz]
        cuty = y[icutz]
        cutz = z[icutz]
        sminmin = 1.e10 * np.ones_like(cutx)
    else:
        # Time to return
        return nea
    # Find closest distance to each arm and then assign arm
    #   We will brute force with a spline
    #   Could get memory intensive
    #   Could refine
    for j in range(adict['narms']):
        # do 50 j=1,narms
        jj = adict['armmap'][j]
        # Crude
        xarm = adict['arm'][j, :adict['kmax'][j], 0]
        yarm = adict['arm'][j, :adict['kmax'][j], 1]

        xtmp = np.outer(cutx, np.ones_like(xarm))
        ytmp = np.outer(cuty, np.ones_like(yarm))
        xatmp = np.outer(np.ones_like(cutx), xarm)
        yatmp = np.outer(np.ones_like(cuty), yarm)
        dist = (xtmp-xatmp)**2 + (ytmp-yatmp)**2
        kmin = np.argmin(dist, axis=1)

        # Refine with a Spline? -- Requires a loop
        csp_xarm = CubicSpline(np.arange(adict['kmax'][j]),
                               adict['arm'][j, :adict['kmax'][j], 0])
        csp_yarm = CubicSpline(np.arange(adict['kmax'][j]),
                               adict['arm'][j, :adict['kmax'][j], 1])
        jtmp = np.zeros((len(cutx), nfine))
        for ii, ix in enumerate(cutx):
            j0 = max(0, kmin[ii]-1)
            j1 = min(adict['kmax'][j], kmin[ii]+1)
            jtmp[ii, :] = np.linspace(j0, j1, num=nfine)
        xjarm = csp_xarm(jtmp)
        yjarm = csp_yarm(jtmp)
        xtmp = np.outer(cutx, np.ones(nfine))
        ytmp = np.outer(cuty, np.ones(nfine))
        dist = (xtmp-xjarm)**2 + (ytmp-yjarm)**2
        # kmin2 = np.argmin(dist, axis=1)
        min_dist2 = np.amin(dist, axis=1)

        smin = np.sqrt(min_dist2)  # Distance of (x,y,z) from this arm's axis
        # Close enough?
        gd_wa = np.where(smin < (3*wa))[0]
        if len(gd_wa) > 0:
            # ga
            ga = np.exp(-(smin[gd_wa]/warms[jj - 1]/wa)**2)  # arm, get the arm weighting factor
            # Galactocentric radial dependence of arms
            tmp_rr = rr[icutz[gd_wa]]
            lg_rr = tmp_rr > Aa
            if np.sum(lg_rr) > 0:
                ga[lg_rr] *= 1. / (np.cosh((tmp_rr[lg_rr]-Aa)/2.0))**2

            # arm3 reweighting:
            if adict['armmap'][j] == 3:
                th3a = 320.
                th3b = 390.
                th3b = 370.
                th3a = 290.
                th3b = 363.
                th3b = 363.
                fac3min = 0.0
                test3 = thxy[icutz[gd_wa]]-th3a
                neg_t3 = test3 < 0
                test3[neg_t3] += 360.
                gd_t3 = (0. <= test3) & (test3 < (th3b-th3a))
                if np.sum(gd_t3) > 0:
                    arg = 6.2831853*(thxy[icutz[gd_wa]][gd_t3]-th3a)/(th3b-th3a)
                    fac = (1.+fac3min + (1.-fac3min)*np.cos(arg))/2.
                    fac = fac**4.0
                    # Update ga
                    ga[gd_t3] *= fac

            # arm2 reweighting:
            if adict['armmap'][j] == 2:
                th2a = 35.
                th2b = 55.
                test2 = thxy[icutz[gd_wa]]-th2a
                fac = 1.
                if False:
                    #    first: as in tc93 (note different definition of theta)
                    gd_t2_first = (0 <= test2) & (test2 < (th2b-th2a))
                    if np.sum(gd_t2_first) > 0:
                        fac = 1. + test2[gd_t2_first]/(th2b-th2a)
                        fac = 1.		# !**** note turned off
                        # ga=ga*fac
                    gd_t2_first2 = test2 > (th2b-th2a)
                    if np.sum(gd_t2_first2) > 0:
                        fac = 2.
                        fac = 1.  # !**** note turned off
                        # ga=ga*fac
                # c    second:  weaken the arm in a short range:
                th2a = 340.
                th2b = 370.
                # c note fix does nothing if fac2min = 1.0
                fac2min = 0.1
                neg_t2 = test2 < 0.
                test2[neg_t2] += 360.
                gd_t2_snd = (0. <= test2) & (test2 < (th2b-th2a))
                if np.sum(gd_t2_snd) > 0:
                    arg = 6.2831853*(thxy[icutz[gd_wa]][gd_t2_snd]-th2a)/(th2b-th2a)
                    fac = (1.+fac2min + (1.-fac2min)*np.cos(arg))/2.
                    # Update ga
                    ga[gd_t2_snd] *= fac

            # Find the ones which are closer than any previous arm
            iclosest = smin[gd_wa] < sminmin[gd_wa]
            if np.sum(iclosest) > 0:
                whicharm[icutz[gd_wa[iclosest]]] = adict['armmap'][j]

            # Update nea
            nea[icutz[gd_wa]] += narms[jj - 1] * (
                ga / (np.cosh(cutz[gd_wa]/(harms[jj - 1]*ha))**2))
    if flg_float:
        nea = float(nea[0])
        whicharm = int(whicharm[0])

    # Return
    return nea

    '''
    Farms = 0
    if(whicharm_spiralmodel .eq. 0) then
    whicharm = 0
    else
      whicharm = armmap(whicharm_spiralmodel)	! remap arm number
      Farms = Fa * farm(whicharm)
    endif
    return
    end
    '''
