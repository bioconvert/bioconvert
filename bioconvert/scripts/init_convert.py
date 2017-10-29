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
""".. rubric:: Create template of bioconvert converter Python file"""
import os
import sys
import argparse


def main(args=None):

    if args is None:
        args = sys.argv[:]

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



