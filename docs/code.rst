====
Code
====

This document describes the code.

NE2001
++++++

Overall, the current implementation is a nearly direct port
of the FORTRAN code created and kindly provide by Cordes & Lazio
(`astro-ph/3101598 <http://adsabs.harvard.edu/abs/2003astro.ph..1598C>`_)

It is important to appreciate, however, that the code is not fully
consistent with the detailed description provided in the publication.
These are the key differences:

* The electron density near the Galactic Center is a constant value
* In NE2001, the thin disk assumes a height of 1.8kpc instead of the A2 parameter
* In NE2001, the Galactic parameters are substantially different than those in the posting

In addition, we have made a few modifications of our own:

* In the original NE2001 code for clumps, q2 <= 5 is used instead of q <= 5.
