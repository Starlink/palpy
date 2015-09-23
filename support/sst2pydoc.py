"""

Read the SST (Starlink Software Tool) prologue(s) from a file and
make the contents available as a dict().

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

import re


def read_prologs(filename):
    """
    Given a filename, search for SST prologues
    and returns a dict where the keys are the name of the
    prolog (from the "Name" field) and the keys are another
    dict with keys correspding fo the SST labels in lowercase.

    prologs = read_prologs( filename )

    Common SST labels are: name, purpose, lanugage, invocation,
    arguments, description, authors, notes etc.
    """
    results = {}
    prolog = {}
    heading_re = re.compile(r"^\*  ([A-Z].*):$")
    heading = ""
    content = ""
    counter = 0

    for line in open(filename):
        line = line.strip()

        # Start of a completely new prolog so reset everything
        if line.startswith("*+"):
            if counter != 0:
                raise ValueError("Started prologue without closing previous prologue")
            prolog = {}
            heading = ""
            content = ""
            counter = counter + 1
            continue

        # End of a prolog. Must store the current dict
        if line.startswith("*-"):
            counter = 0
            if len(heading):
                # Flush current heading
                prolog[heading] = content
                content = ""
            name = prolog['name'].strip()
            results[name] = prolog
            prolog = None
            continue

        # If we are not in a prologue then nothing further is needed
        if counter == 0:
            continue

        counter = counter + 1

        # Completely blank lines are ignored
        if len(line) == 0:
            continue

        # Look for a new section heading
        match_head = heading_re.search(line)
        if match_head is not None:
            if len(heading):
                # Flush previous heading
                prolog[heading] = content
            heading = match_head.group(1).lower()
            content = ""
            continue

        if line.startswith("*     "):
            content = content + line[6:] + "\n"
            continue
        elif line == "*":
            content = content + "\n"
            continue

        if counter:
            raise ValueError("Error parsing SST prologue line "+str(counter)+":'" + line + "'")
    return results

if __name__ == "__main__":
    results = read_prologs("cextern/pal/palAddet.c")
    print(results)
