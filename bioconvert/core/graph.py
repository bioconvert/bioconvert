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
import os
import json
import glob
from os.path import join, basename


def create_graph(filename, layout="dot", use_singularity=False):
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
            raise
        from pygraphviz import AGraph
        dg = AGraph(directed=True)

        for a, b in rr.get_conversions():
            dg.add_edge(a, b)

        dg.layout(layout)
        dg.draw(filename)

    except:

        dot = """
strict digraph{
    node [label="\\N"];

    """
        nodes =  set([item for items in rr.get_conversions() for item in items])

        for node in nodes:
            dot += "{};\n".format(node)
        for a, b in rr.get_conversions():
            dot+="{} -> {};\n".format(a, b)
        dot += "}\n"

        from easydev import TempFile
        from bioconvert import shell
        dotfile = TempFile(suffix=".dot")
        with open(dotfile.name, "w") as fout:
            fout.write(dot)

        dotpath = ""
        if use_singularity:
            # download singularity
            from bioconvert import configuration as config
            from easydev import md5
            # note that in singularity v2.4, whatever extension you put, it is
            # replaced by simg
            singfile = "{}/graphviz.simg".format(config.user_config_dir)

            print(singfile)

            if os.path.exists(singfile) and md5(singfile) == "4288088d91c848e5e3a327282a1ab3d1":
                print("Found singularity (graphviz) image")
            else:
                print("Downloading singularity. Please wait")
                cmd = "singularity pull --name {}  shub://cokelaer/graphviz4all:v1"
                cmd = cmd.format(singfile)
                try:
                    shell(cmd)
                except:
                    import os
                    os.system(cmd)

            dotpath = "singularity run {} ".format(singfile)

        ext = filename.rsplit(".", 1)[1]
        cmd = "{}dot -T{} {} -o {}".format(dotpath, ext, dotfile.name, filename)
        try:
            shell(cmd)
        except:
            import os
            os.system(cmd)
        dotfile.delete()
















