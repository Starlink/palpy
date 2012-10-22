PALpy - Python positional astronomy library
===========================================

PALpy is a python interface to the PAL positional astronomy library
(<https://github.com/Starlink/pal>). All source code is included with
the distribution via git submodules.

To build the python interface from the git repository:

    git submodule init
    git submodule update
    python setup.py build
    python setup.py install
    python test_pal.py

Requires that cython and numpy are installed.

The interface uses the same names as the C interface with the pal
prefix missing:

    import palpy as pal

    djm = pal.caldj( 1999, 12, 31 )
    (r2000, d2000) = pal.fk45z( 1.2, -0.3, 1960 )

All arguments modified by the C API are returned. No arguments
are modified. The one routine that is different is palObs which
returns a simple dict that can be searched using standard python.
The keys to the dict are the short names and the values are another
dict with keys name, long, lat and height.

See `test_pal.py` for detailed examples of all functions.

Further Work
------------

Currently a handful of PAL routines are missing from the python interface.
They can trivially be added when time is available and there is a demand.
Furthermore not all the SLALIB routines are available from PAL.

Author
------

Tim Jenness (tim.jenness@gmail.com).

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 2 of
the License, or (at your option) any later version.

This program is distributed in the hope that it will be
useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA.
