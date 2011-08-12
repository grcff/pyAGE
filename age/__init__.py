# age/__init__.py
#  pyAGE - A Python implementation of the Analog Genetic Encoding
#  Copyright (C) 2010  Janosch Gräf
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.


def __load_doc__(filename = "__doc__.txt"):
    import os
    global __path__

    path = os.path.join(__path__[0], filename)
    try:
        f = open(path, "rt")
        doc = f.read()
        f.close()
        return doc
    except IOError:
        return "Could not load documentation: "+path+"\n"

__version__ = "0.1"
__doc__ = __load_doc__()
__all__ = ["Descriptor", "Genome", "Device", "Population", "Agent"]

# import classes
from .descriptor import *
from .genome import *
from .population import *
