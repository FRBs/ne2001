=======
Scripts
=======

This document describes the scripts packaged
with ne2001.

ne2001_lb
+++++++++

Calculate quantities along an input Galactic
line-of-sight.  Main inputs are longitude and
latitude.  Here is the usage::

    ne2001_lb -h
    usage: ne2001_lb [-h] [-d D] l b

    Calculate quantities along a Galactic sightline v0.1

    positional arguments:
      l           Galactic longitude (deg)
      b           Galactic latitude (deg)

    optional arguments:
      -h, --help  show this help message and exit
      -d D        Distance (kpc)

Currently, this script calculates the DM along the sightline
to a default distance of 100 kpc.
