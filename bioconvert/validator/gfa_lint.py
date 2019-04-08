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
import colorlog
_log = colorlog.getLogger(__name__)


__all__ = ['GFALint']


class GFALint(object):
    """

    see https://github.com/sjackman/gfalint/
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
                if line[0] in "#EFGHLOPSU":
                    pass
                elif len(line.strip()) == 0:
                    _log.warning("Found empty line on line {}".format(line))
                else:
                    raise ValueError("Unknown starting field ({}) on line {}".format((line[0], i)))

