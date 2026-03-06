import pytest

from bioconvert.pajek2graphml import PAJEK2GRAPHML
from bioconvert import TempFile

from . import test_dir


@pytest.mark.parametrize("method", PAJEK2GRAPHML.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/pajek/test_network.net"
    with TempFile(suffix=".graphml") as outfile:
        convert = PAJEK2GRAPHML(infile, outfile.name)
        convert(method=method)

        import networkx as nx

        G_in = nx.read_pajek(infile)
        G_out = nx.read_graphml(outfile.name)
        assert G_in.number_of_nodes() == G_out.number_of_nodes()
        assert G_in.number_of_edges() == G_out.number_of_edges()
