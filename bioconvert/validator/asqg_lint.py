# -*- coding: utf-8 -*-
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
"""ASQG validator"""
import colorlog
_log = colorlog.getLogger(__name__)


class ASQGLint():
    """ASQG format validator

    """
    def __init__(self, filename):
        self.filename = filename


    def validate(self):
        # read line by line. Checks 
        # - lines start with HT, VT or ED
        # - lines must be tab delimited
        # - lines VT must have 2 fields only
        with open(self.filename, "r") as fh:
            for i, line in enumerate(fh.readlines()):
                if line[0:2] == "HT":
                    pass
                elif line[0:2] == "VT":
                    if len(line.split()) < 3:
                        raise ValueError("Incorrect format on line {}. ".format(i) +
                            "Line starting with VT must have at least 3 fields")
                elif line[0:2] == "ED":
                    if len(line.split()) != 11:
                        raise ValueError("Incorrect format on line {}. ".format(i) +
                            "Line starting with ED must have 10 fields + ED prefix")
                elif len(line.strip()) == 0:
                    _log.warning("Found empty line on line {}".format(line))
                else:
                    raise ValueError()

