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
""".. rubric:: Create template of bioconvert converter Python file"""
import os
import sys
import argparse


def main(args=None):

    if args is None:
        args = sys.argv[1:]


    from easydev.console import purple, underline

    arg_parser = argparse.ArgumentParser(prog="bioconvert_init",
                                         epilog=" ----    ",
                                         description="""DESCRIPTION:

Create a Python module to ease addition of new converters

""")
    arg_parser.add_argument("-i", "--input-extension",
                            help="input_extension")
    arg_parser.add_argument("-o", "--output-extension",
                            help="output_extension")

    args = arg_parser.parse_args()


    from bioconvert.core.init import InitConverter
    ic = InitConverter(args.input_extension, args.output_extension)
    print(ic.get_content())


if __name__ == "__main__":
    main()



