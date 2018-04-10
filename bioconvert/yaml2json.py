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
"""Convert :term:`YAML` to :term:`JSON` format"""
import yaml, json
from bioconvert import ConvBase
import colorlog

from bioconvert.core.decorators import requires_nothing

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

    def get_json(self):
        """Return the JSON dictionary corresponding to the YAML input"""
        data = yaml.load(open(self.infile, "r"))
        return json.dumps(data, sort_keys=True, indent=4)

    @requires_nothing
    def _method_python(self, *args, **kwargs):
        with open(self.outfile, "w") as outfile:
            outfile.write(self.get_json())
