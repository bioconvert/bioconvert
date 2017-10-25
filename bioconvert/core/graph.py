import json
import glob
from os.path import join, basename
# install this with "conda install -c conda-forge python-graphviz"
import pygraphviz as pgv



def create_graph(filename, layout="dot"):
    """

    :param filename: should end in .png or .svg
    """
    from bioconvert.core.registry import Registry
    rr = Registry()

    dg = pgv.AGraph()

    for a, b in rr.get_conversions():
        dg.add_edge(a, b)


    dg = dg.to_directed()
    dg.layout(layout)
    dg.draw(filename)

