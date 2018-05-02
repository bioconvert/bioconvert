from bioconvert.core.graph import create_graph, create_graph_for_cytoscape
from easydev import TempFile


def test_create_graph_singularity():
    with TempFile(suffix=".png") as fout:
        create_graph(fout.name, use_singularity=True)


def test_create_graph():
    with TempFile(suffix=".png") as fout:
        create_graph(fout.name, use_singularity=False)


def test_create_cytoscape_export():
    create_graph_for_cytoscape(True)
    create_graph_for_cytoscape(False)
