# -*- coding: utf-8 -*-

###########################################################################
# Bioconvert is a project to facilitate the interconversion               #
# of life science data from one format to another.                        #
#                                                                         #
# Copyright © 2018-2022  Institut Pasteur, Paris and CNRS.                     #
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

"""
Available methods per converter
=====================================

Plot number of implemented methods per converter.


"""
#################################################
#
from bioconvert.core.registry import Registry

r = Registry()
info = r.get_info()

# The available unique converters
converters = [x for x in info.items()]

# the number of methods per converter
data = [info[k] for k, v in info.items()]

# the number of formats
A1 = [y for x in list(r.get_conversions()) for y in x[0]]
A2 = [y for x in list(r.get_conversions()) for y in x[1]]
formats = set(A1 + A2)

print("Number of formats: {}".format(len(formats)))
print("Number of converters: {}".format(len(converters)))
print("Number of methods : {}".format(sum(data)))

#####################################################
from pylab import hist, clf, xlabel, grid

clf()
hist(data, range(17), ec="k", zorder=2, align="left")
xlabel("Number of methods")
grid(zorder=-1)
