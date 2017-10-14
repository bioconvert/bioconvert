import json
import glob
from os.path import join, basename
# install this with "conda install -c conda-forge python-graphviz"
import graphviz as gv



def create_graph(filename, format="svg"):
    from bioconvert.core.registry import Registry
    rr = Registry()
    dg = gv.Digraph(filename=filename, format=format)

    for a,b in rr.get_conversions():
        dg.node(a)
        dg.node(b)
        dg.edge(a, b)
    dg.render()
