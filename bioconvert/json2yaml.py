# -*- coding: utf-8 -*-
#
#  This file is part of Bioconvert software
#
#  Copyright (c) 2017 - Bioconvert Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/biokit/bioconvert
#  documentation: http://bioconvert.readthedocs.io
#
##############################################################################
"""Convert :term:`JSON` to :term:`YAML` format"""
import yaml, json
from bioconvert import ConvBase
from bioconvert.core.decorators import requires, requires_nothing

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
    def _method_yaml(self, *args, **kwargs):
        with open(self.infile, "r") as infile:
            data = json.load(infile)
        with open(self.outfile, "w") as outfile:
            yamldata = yaml.dump(data, Dumper=yaml.dumper.Dumper,
                default_flow_style="", indent=4)
            outfile.write(yamldata)

