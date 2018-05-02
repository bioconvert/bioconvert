from bioconvert.core.graph import create_graph, create_graph_for_cytoscape
from easydev import TempFile


def test_create_graph_singularity():
    with TempFile(suffix=".png") as fout:
        create_graph(fout.name, use_singularity=True)


def test_create_graph():
    with TempFile(suffix=".png") as fout:
        create_graph(fout.name, use_singularity=False)


def test_create_cytoscape_export():
    g_small = create_graph_for_cytoscape(False)
    assert len(g_small["elements"]["nodes"]) > 0
    assert len(g_small["elements"]["edges"]) > 0
    g_full = create_graph_for_cytoscape(True)
    assert len(g_full["elements"]["nodes"]) >= len(g_small["elements"]["nodes"])
    assert len(g_full["elements"]["edges"]) >= len(g_small["elements"]["edges"])
