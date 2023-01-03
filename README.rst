PALpy - Python positional astronomy library
===========================================

.. image:: https://zenodo.org/badge/12517/Starlink/palpy.svg
   :target: https://zenodo.org/badge/latestdoi/12517/Starlink/palpy

.. image:: https://travis-ci.org/Starlink/palpy.svg
    :target: https://travis-ci.org/Starlink/palpy

PALpy is a python interface to the PAL positional astronomy library
(<https://github.com/Starlink/pal>). All source code is included with
the distribution via git submodules.

To build the python interface from the git repository::

    git submodule init
    git submodule update
    python setup.py build
    python setup.py install
    python test_pal.py

Requires that Cython and numpy are installed.

The interface uses the same names as the C interface with the pal
prefix missing::

    import palpy as pal

    djm = pal.caldj( 1999, 12, 31 )
    (r2000, d2000) = pal.fk45z( 1.2, -0.3, 1960 )

All arguments modified by the C API are returned. No arguments
are modified. The one routine that is different is palObs which
returns a simple dict that can be searched using standard python.
The keys to the dict are the short names and the values are another
dict with keys name, long, lat and height.

Docstrings are included with each routine describing the python
usage and the associated C documentation:

     help(palpy)
     help(palpy.caldj)

See ``test_pal.py`` for detailed examples of all functions.

Documentation:
--------------

The module provides documentation on how to use the python interface
to PAL. It does not contain information on how to use
PAL itself.

The basic PAL documentation can be found at

 http://www.starlink.ac.uk/star/docs/sun267.htx/sun267.html

For more information the SLA documentation provides much more
detail:

  http://www.starlink.ac.uk/star/docs/sun67.htx/sun67.html

A description paper for PAL was published in the ADASS XXII
conference proceedings:

  http://adsabs.harvard.edu/abs/2013ASPC..475..307J

Please consider citing this if you use PAL or palpy in your
research.

Further Work
------------

Not all the SLALIB routines are available from PAL. More can be added
as need arises.

Author
------

Copyright (C) 2012, 2014-2016
Tim Jenness (tim.jenness@gmail.com).
All Rights Reserved.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 3 of
the License, or (at your option) any later version.

This program is distributed in the hope that it will be
useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
