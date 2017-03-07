=====
Usage
=====

Overview
++++++++

To use ne2001 in a project::

	import ne2001



Simple Calculations
+++++++++++++++++++

Here is an example for calculating the DM along a Galactic
sightline to a distance of 20kpc::

    from ne2001 import io as ne_io
    PARAMS = ne_io.read_params()
    ne = density.ElectronDensity(**PARAMS)
    l, b  = 0., 90.
    DM = ne.DM(l, b, 20.)
