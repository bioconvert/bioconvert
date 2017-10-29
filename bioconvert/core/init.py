# -*- coding: utf-8 -*-
##############################################################################
#  This file is part of Bioconvert software
#
#  Copyright (c) 2017 - Bioconvert Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/biokit/bioconvert
#  documentation: http://bioconvert.readthedocs.io
##############################################################################
import os



template = '''

# -*- coding: utf-8 -*-
##############################################################################
#  This file is part of Bioconvert software
#
#  Copyright (c) 2017 - Bioconvert Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/biokit/bioconvert
#  documentation: http://bioconvert.readthedocs.io
##############################################################################
""" description """

from bioconvert import ConvBase

__all__ = ["{input}2{output}"]


class {input}2{output}(ConvBase):
    """Convert :term:`{input}` file to :term:`{output}` file

    Some description.

    """
    input_ext = [".{inputext}"]
    output_ext = [".{outputext}"]

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input {input} file
        :param str outfile: output {output} filename

        """
        super({input}2{output}, self).__init__(infile, outfile, *args, **kargs)

        self._default_method = "default"

    def _method_default(self, *args, **kwargs):
        """some description"""
        cmd = "ls "
        # use self.infile, self.outfile
        self.execute(cmd)
'''


class InitConverter():
    def __init__(self, inputext, outputext):

        self.input = inputext
        self.output = outputext

    def get_content(self):
        return template.format(
                input=self.input.upper(),
                output=self.output.upper(),
                inputext=self.input,
                outputext=self.output)












