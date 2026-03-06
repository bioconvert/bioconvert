import pytest

from bioconvert.graphml2pajek import GRAPHML2PAJEK
from bioconvert import TempFile

from . import test_dir


@pytest.mark.parametrize("method", GRAPHML2PAJEK.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/graphml/test_network.graphml"
    with TempFile(suffix=".net") as outfile:
        convert = GRAPHML2PAJEK(infile, outfile.name)
        convert(method=method)

        import networkx as nx

        G_in = nx.read_graphml(infile)
        G_out = nx.read_pajek(outfile.name)
        assert G_in.number_of_nodes() == G_out.number_of_nodes()
        assert G_in.number_of_edges() == G_out.number_of_edges()
