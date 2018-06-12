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

"""
Possible Conversion
====================

Plot directed graph of possible conversion


"""
#################################################
#
from bioconvert.core.graph import create_graph



#####################################################
# If you use pygraphviz, you can have a good quality 
# image using:
import matplotlib as mpl
mpl.rcParams['figure.dpi']= 250

############################################################
# In order to create the following image, you need graphviz
# and pygraphviz.
# If you cannot install those packages, you may use a singularity
# image like in the following example by setting the use_singularity
# parameter to True. This would work under Linux. Not tested on other systems
# yet.
create_graph("conversion.png", use_singularity=True)

#####################################################
#
from pylab import imshow, imread, xticks, yticks, gca
imshow(imread("conversion.png"), interpolation="nearest")
xticks([])
yticks([])
ax = gca()
ax.axis('off')
