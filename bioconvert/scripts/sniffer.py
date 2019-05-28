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
""".. rubric:: Standalone application dedicated to conversion"""
import argparse
import json
import sys

import bioconvert
from bioconvert import ConvBase
from bioconvert.core import graph
from bioconvert.core.base import ConvMeta
from bioconvert.core.converter import Bioconvert
from bioconvert.core.decorators import get_known_dependencies_with_availability
from bioconvert.core.registry import Registry


def main(args=None):
    if args is None:
        args = sys.argv[1:]

    arg_parser = argparse.ArgumentParser(prog="bioconvert_sniffer",
                                         description="""Infer input format (in progress)""",
                                         formatter_class=argparse.RawDescriptionHelpFormatter,
                                         epilog="""
Please visit http://bioconvert.readthedocs.org for more information about the
project or formats available.

Bioconvert is an open source collaborative project. Please feel free to 
join us at https://github/biokit/bioconvert
""")


    arg_parser.add_argument("-v", "--verbosity",
                            default=bioconvert.logger.level,
                            help="Set the outpout verbosity.",
                            choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                            )
    arg_parser.add_argument("--input", default="input",
                            help="Set the input file.")

    args = arg_parser.parse_args(args)


    from bioconvert.io.sniffer import Sniffer
    print(args)

    s = Sniffer()
    ret = s.sniff(args.input)
    url = "https://bioconvert.readthedocs.io/en/master/formats.html"
    if ret:
        print("\nFound possible candidate(s): {}".format(ret))
        if isinstance(ret, list):
            for this in ret:
                print("{}#{}".format(url, this))
        else:
            print("{}#{}".format(url, ret))
    else:
        print("No candidate format identified. This is a tool in progress as"
    "shown in the github entry https://github.com/bioconvert/bioconvert/issues/225")



if __name__ == "__main__":
    main()
