import pytest

from bioconvert.graphml2gml import GRAPHML2GML
from bioconvert import TempFile

from . import test_dir


@pytest.mark.parametrize("method", GRAPHML2GML.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/graphml/test_network.graphml"
    with TempFile(suffix=".gml") as outfile:
        convert = GRAPHML2GML(infile, outfile.name)
        convert(method=method)

        import networkx as nx

        G_in = nx.read_graphml(infile)
        G_out = nx.read_gml(outfile.name)
        assert G_in.number_of_nodes() == G_out.number_of_nodes()
        assert G_in.number_of_edges() == G_out.number_of_edges()
