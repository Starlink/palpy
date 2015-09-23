"""

Extract the upstream PAL version number from the PAL
C source distribution.

Authors
-------

Tim Jenness

Licence
-------

This program is free software: you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General
License along with this program.  If not, see
<http://www.gnu.org/licenses/>.

"""

import os
import re


def read_pal_version():
    """
    Scans the PAL configure.ac looking for the version number.

    (vers, maj, min, patchlevel) = read_pal_version()

    Returns the version as a string and the major, minor
    and patchlevel version integers

    """
    verfile = os.path.join("cextern", "pal", "configure.ac")
    verstring = "-1.-1.-1"
    for line in open(verfile):
        if line.startswith("AC_INIT"):
            # Version will be in string [nn.mm.pp]
            match = re.search(r"\[(\d+\.\d+\.\d+)\]", line)
            if match:
                verstring = match.group(1)
            break
    (major, minor, patch) = verstring.split(".")
    return (verstring, major, minor, patch)

if __name__ == "__main__":
    v, maj, min, p = read_pal_version()
    print(v, maj, min, p)
