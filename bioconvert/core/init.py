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
from bioconvert import requires

__all__ = ["{input}2{output}"]


class {input}2{output}(ConvBase):
    """Convert :term:`{input}` file to :term:`{output}` file

    Some description.

    """

    _default_method = "default"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input {input} file
        :param str outfile: output {output} filename

        """
        super({input}2{output}, self).__init__(infile, outfile, *args, **kargs)

    @requires("ls")
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












