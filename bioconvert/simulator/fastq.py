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


class FastqSim():
    def __init__(self, outfile):
        self.outfile = outfile
        self.nreads = 1000000
        # do we want a wrap or unwrap version ?
        # For now, unwrap
        self.wrap = False
        self.read_length = 250

    def simulate(self):
        RL = self.read_length
        with open(self.outfile, "w") as fout:
            for i in range(self.nreads):
                fout.write("@identifier whatever it means but long enough\n")
                fout.write("ACGT" * (RL // 4) + "A" * (RL % 4) + "\n")
                fout.write("+\n")
                fout.write("C" * RL + "\n")

