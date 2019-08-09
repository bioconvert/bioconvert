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
"""Convert :term:`JSON` to :term:`YAML` format"""
import yaml
import json
from bioconvert import ConvBase
from bioconvert.core.decorators import requires, requires_nothing, compressor


__all__ = ["JSON2YAML"]


class JSON2YAML(ConvBase):
    """Convert :term:`JSON` file into :term:`YAML` file

    Conversion is based on yaml and json standard Python modules
    Indentation is set to 4 by default and affects the sections (not the list).
    For example::

        fruits_list:
        - apple
        - orange
        section1:
            do: true
            misc: 1

    """

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input JSON file
        :param str outfile: input YAML file. 
        """
        super(JSON2YAML, self).__init__(infile, outfile, *args, **kargs)

    @requires_nothing
    @compressor
    def _method_yaml(self, *args, **kwargs):
        with open(self.infile, "r") as infile:
            data = json.load(infile)
        with open(self.outfile, "w") as outfile:
            yamldata = yaml.dump(data, Dumper=yaml.dumper.Dumper,
                default_flow_style="", indent=4)
            outfile.write(yamldata)

