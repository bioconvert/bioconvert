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

    arg_parser = argparse.ArgumentParser(prog="bioconvert",
                                         description="""Convertor infer the
                                         formats from the first command. We do
                                         not scan the input file. Therefore
                                         users must ensure that their input
                                         format files are properly
                                         formatted.""",
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
    arg_parser.add_argument("--no-plot", action="store_true")

    args = arg_parser.parse_args(args)


    from bioconvert.core.registry import Registry
    r = Registry()
    info = r.get_info()

    # The available unique converters
    converters = [x for x in info.items()]

    # the number of methods per converter
    data = [info[k] for k,v in info.items()]

    # the number of formats
    A1 = [x[0] for x in list(r.get_conversions())]
    A2 = [x[1] for x in list(r.get_conversions())]
    formats = set(A1 + A2)


    print("Number of formats: {}".format(len(formats)))
    print("Number of converters: {}".format(len(converters)))
    print("Number of methods : {}".format(sum(data)))

    if args.no_plot is False:
        from pylab import hist, clf, xlabel, grid, show
        clf()
        hist(data, range(17), ec="k", zorder=2, align="left")
        xlabel("Number of methods")
        grid(zorder=-1)
        show()



if __name__ == "__main__":
    main()
