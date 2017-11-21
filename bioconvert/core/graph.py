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
""" description """
import json
import glob
from os.path import join, basename


def create_graph(filename, layout="dot"):
    """

    :param filename: should end in .png or .svg
    """
    from bioconvert.core.registry import Registry
    rr = Registry()

    try:
        from pygraphviz import AGraph
    except:
        print("Please install graphviz/pygraphviz")
        class AGraph():
            pass

    dg = AGraph(directed=True)

    for a, b in rr.get_conversions():
        dg.add_edge(a, b)

    dg.layout(layout)
    dg.draw(filename)

