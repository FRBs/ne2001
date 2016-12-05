========
Overview
========

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - | |travis| |requires|
        | |coveralls| |codecov|
        | |landscape| |scrutinizer| |codacy| |codeclimate|
..    * - package
..      - |version| |downloads| |wheel| |supported-versions| |supported-implementations|

   
.. |docs| image:: https://readthedocs.org/projects/ne2001/badge/?style=flat
    :target: https://readthedocs.org/projects/ne2001
    :alt: Documentation Status

.. |travis| image:: https://travis-ci.org/benbaror/ne2001.svg?branch=master
    :alt: Travis-CI Build Status
    :target: https://travis-ci.org/benbaror/ne2001

.. |requires| image:: https://requires.io/github/benbaror/ne2001/requirements.svg?branch=master
    :alt: Requirements Status
    :target: https://requires.io/github/benbaror/ne2001/requirements/?branch=master

.. |coveralls| image:: https://coveralls.io/repos/benbaror/ne2001/badge.svg?branch=master&service=github
    :alt: Coverage Status
    :target: https://coveralls.io/r/benbaror/ne2001

.. |codecov| image:: https://codecov.io/github/benbaror/ne2001/coverage.svg?branch=master
    :alt: Coverage Status
    :target: https://codecov.io/github/benbaror/ne2001

.. |landscape| image:: https://landscape.io/github/benbaror/ne2001/master/landscape.svg?style=flat
    :target: https://landscape.io/github/benbaror/ne2001/master
    :alt: Code Quality Status

.. |codacy| image:: https://img.shields.io/codacy/REPLACE_WITH_PROJECT_ID.svg?style=flat
    :target: https://www.codacy.com/app/benbaror/ne2001
    :alt: Codacy Code Quality Status

.. |codeclimate| image:: https://codeclimate.com/github/benbaror/ne2001/badges/gpa.svg
   :target: https://codeclimate.com/github/benbaror/ne2001
   :alt: CodeClimate Quality Status

.. |version| image:: https://img.shields.io/pypi/v/ne2001.svg?style=flat
    :alt: PyPI Package latest release
    :target: https://pypi.python.org/pypi/ne2001

.. |downloads| image:: https://img.shields.io/pypi/dm/ne2001.svg?style=flat
    :alt: PyPI Package monthly downloads
    :target: https://pypi.python.org/pypi/ne2001

.. |wheel| image:: https://img.shields.io/pypi/wheel/ne2001.svg?style=flat
    :alt: PyPI Wheel
    :target: https://pypi.python.org/pypi/ne2001

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/ne2001.svg?style=flat
    :alt: Supported versions
    :target: https://pypi.python.org/pypi/ne2001

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/ne2001.svg?style=flat
    :alt: Supported implementations
    :target: https://pypi.python.org/pypi/ne2001

.. |scrutinizer| image:: https://img.shields.io/scrutinizer/g/benbaror/ne2001/master.svg?style=flat
    :alt: Scrutinizer Status
    :target: https://scrutinizer-ci.com/g/benbaror/ne2001/


   
.. end-badges

Python implementation of Cordes-Lazio's NE2001 Galactic Free Electron Density Model

* Free software: BSD license

.. Installation
.. ============

.. ::

..    pip install ne2001

.. Documentation
.. =============

.. https://ne2001.readthedocs.io/

Development
===========

To run the all tests run::

    tox

Note, to combine the coverage data from all the tox environments run:

.. list-table::
    :widths: 10 90
    :stub-columns: 1

    - - Windows
      - ::

            set PYTEST_ADDOPTS=--cov-append
            tox

    - - Other
      - ::

            PYTEST_ADDOPTS=--cov-append tox
