# -*- coding: utf-8 -*-
#!/opt/local/bin/python

###########################################################################
# Bioconvert is a project to facilitate the interconversion               #
# of life science data from one format to another.                        #
#                                                                         #
# Authors: see CONTRIBUTORS.rst                                           #
# Copyright © 2018  Institut Pasteur, Paris and CNRS.                     #
# See the COPYRIGHT file for details                                      #
#                                                                         #
# bioconvert is free software: you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or       #
# (at your option) any later version.                                     #
#                                                                         #
# bioconvert is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of          #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
# GNU General Public License for more details.                            #
#                                                                         #
# You should have received a copy of the GNU General Public License       #
# along with this program (COPYING file).                                 #
# If not, see <http://www.gnu.org/licenses/>.                             #
###########################################################################

import sys
import mmap

with open(sys.argv[1], "r+") as inp:
    with open(sys.argv[2], "wb") as out:
        mapp = mmap.mmap(inp.fileno(), 0)
        line = mapp.readline()
        while line:
            out.write(b">")
            out.write(line[1:])
            out.write(mapp.readline())
            mapp.readline()
            mapp.readline()
            line = mapp.readline()
        mapp.close()