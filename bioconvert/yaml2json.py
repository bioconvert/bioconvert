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
"""Convert :term:`YAML` to :term:`JSON` format"""
import yaml
import json
from bioconvert import ConvBase
import colorlog

from bioconvert.core.decorators import requires_nothing
from bioconvert.core.decorators import compressor

logger = colorlog.getLogger(__name__)

__all__ = ["YAML2JSON"]


class YAML2JSON(ConvBase):
    """Convert :term:`YAML` file into :term:`JSON` file

    Conversion is based on yaml and json standard Python modules

    .. note:: YAML comments will be lost in JSON output

    :reference: http://yaml.org/spec/1.2/spec.html#id2759572
    """

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input YAML file.
        :param str outfile: input JSON file
        """
        super(YAML2JSON, self).__init__(infile, outfile, *args, **kargs)

    @requires_nothing
    @compressor
    def get_json(self):
        """Return the JSON dictionary corresponding to the YAML input"""
        data = yaml.load(open(self.infile, "r"))
        return json.dumps(data, sort_keys=True, indent=4)

    @requires_nothing
    @compressor
    def _method_python(self, *args, **kwargs):
        with open(self.outfile, "w") as outfile:
            outfile.write(self.get_json())