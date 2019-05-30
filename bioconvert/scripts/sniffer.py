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
import os

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

USAGE:

    bioconvert_sniffer test.bam
    bioconvert_sniffer --input test.bam --verbosity INFO

""")


    arg_parser.add_argument("-v", "--verbosity",
                            default=bioconvert.logger.level,
                            help="Set the outpout verbosity.",
                            choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                            )
    arg_parser.add_argument("--input", default="input", nargs="+",
                            help="Set the input file.")
    arg_parser.add_argument("--quiet", default="quiet", action="store_true",
                            help="simply return a two columns output with filename and found format.")


    if args[0][0] != "-":
        args = ["--input"] +  args

    args = arg_parser.parse_args(args)

    from bioconvert.io.sniffer import Sniffer

    # Set the verbosity
    from bioconvert import logger
    logger.level = args.verbosity


    # Sniff files

    s = Sniffer()

    for filename in args.input:

        ret = s.sniff(filename)
        url = "https://bioconvert.readthedocs.io/en/master/formats.html"
        if ret:
            fname = os.path.split(filename)[1]
            if args.quiet:
                pass
            else:
                print("\nFound possible candidate(s) for {}: {}".format(fname, ret))

            if args.quiet:
                print("{}\t {}".format(fname, ret))
            else:
                if isinstance(ret, list):
                    for this in ret:
                        print("See {}#{}".format(url, this))
                else:
                    print("See {}#{}".format(url, ret))
        else:
            
            fname = os.path.split(filename)[1]
            msg = "No candidate format identified for {} This is a tool in progress as"
            msg += "shown in the github entry https://github.com/bioconvert/bioconvert/issues/225"
            if args.quiet:
                print("{}\t NA".format(fname))
            else:
                print(msg.format(fname))



if __name__ == "__main__":
    main()
