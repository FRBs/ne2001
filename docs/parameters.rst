==========
Parameters
==========

Components of the Galactic ISM are parameterized
in terms of size, central electron density, etc.

Currently, the code reads the parameters from the JSON
file data/ne2001_params.json.  You can modify this
directly and future code will allow easier manipulation.

v0.1
++++

Here is the JSON file::

    {
        "version": "v0.1 by JXP on 05March2017",
        "galactic_center": {
            "F": 60000.0,
            "center": [
                -0.01,
                0.0,
                -0.02
            ],
            "e_density": 10.0,
            "height": 0.026,
            "radius": 0.145
        },
        "ldr": {
            "F": 0.1,
            "center": [
                1.36,
                8.06,
                0.0
            ],
            "e_density": 0.012,
            "ellipsoid": [
                1.5,
                0.75,
                0.5
            ],
            "theta": -0.4223696789826278
        },
        "lhb": {
            "F": 0.01,
            "center": [
                0.01,
                8.45,
                0.17
            ],
            "cylinder": [
                0.085,
                0.1,
                0.33
            ],
            "e_density": 0.005,
            "theta": 0.2617993877991494
        },
        "loop_in": {
            "F": 0.2,
            "center": [
                -0.045,
                8.4,
                0.07
            ],
            "e_density": 0.0125,
            "radius": 0.12
        },
        "loop_out": {
            "F": 0.01,
            "center": [
                -0.045,
                8.4,
                0.07
            ],
            "e_density": 0.0125,
            "radius": 0.18
        },
        "lsb": {
            "F": 0.01,
            "center": [
                -0.75,
                9.0,
                -0.05
            ],
            "e_density": 0.016,
            "ellipsoid": [
                1.05,
                0.425,
                0.325
            ],
            "theta": 2.426007660272118
        },
        "spiral_arms": {
            "e_density": 1.0,
            "F": 1.0
        },
        "thick_disk": {
            "F": 0.18,
            "e_density": 0.034020618556701035,
            "height": 0.97,
            "radius": 17.5
        },
        "thin_disk": {
            "F": 120,
            "e_density": 0.08,
            "height": 0.15,
            "radius": 3.8
        }
    }
