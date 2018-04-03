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
from os import environ

import colorlog

_log = colorlog.getLogger(__name__)

def create_graph(filename, layout="dot", use_singularity=False, color_for_disabled_converter='red'):
    """

    :param filename: should end in .png or .svg or .dot

    If extension is .dot, only the dot file is created.
    This is useful if you have issues installing graphviz.
    If so, under Linux you could use our singularity container
    see github.com/cokelaer/graphviz4all

    """
    from bioconvert.core.registry import Registry
    rr = Registry()

    try:
        if filename.endswith(".dot") or use_singularity is True:
            raise Exception()
        from pygraphviz import AGraph
        dg = AGraph(directed=True)

        for a, b, s in rr.get_all_conversions():
            dg.add_edge(a, b, color='black' if s else color_for_disabled_converter)

        dg.layout(layout)
        dg.draw(filename)

    except Exception as e:
        _log.error(e)
        dot = """
strict digraph{
    node [label="\\N"];

    """
        nodes = set([item for items in rr.get_all_conversions() for item in items[0:1]])

        for node in nodes:
            dot += "\"{}\";\n".format(node)
        for a, b, s in rr.get_all_conversions():
            dot += "\"{}\" -> \"{}\";\n".format(a, b)
        dot += "}\n"

        from easydev import TempFile
        from bioconvert import shell
        dotfile = TempFile(suffix=".dot")
        with open(dotfile.name, "w") as fout:
            fout.write(dot)

        dotpath = ""
        if use_singularity:
            from bioconvert.core.downloader import download_singularity_image
            singfile = download_singularity_image(
                "graphviz.simg",
                "shub://cokelaer/graphviz4all:v1",
                "4288088d91c848e5e3a327282a1ab3d1")

            dotpath = "singularity run {} ".format(singfile)
            on_rtd = environ.get('READTHEDOCS', None) == 'True'
            if on_rtd:
                dotpath = ""


        ext = filename.rsplit(".", 1)[1]
        cmd = "{}dot -T{} {} -o {}".format(dotpath, ext, dotfile.name, filename)
        try:
            shell(cmd)
        except:
            import os
            os.system(cmd)
        #dotfile.delete()
