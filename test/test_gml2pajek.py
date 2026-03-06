import pytest

from bioconvert.gml2pajek import GML2PAJEK
from bioconvert import TempFile

from . import test_dir


@pytest.mark.parametrize("method", GML2PAJEK.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/gml/test_network.gml"
    with TempFile(suffix=".net") as outfile:
        convert = GML2PAJEK(infile, outfile.name)
        convert(method=method)

        import networkx as nx

        G_in = nx.read_gml(infile)
        G_out = nx.read_pajek(outfile.name)
        assert G_in.number_of_nodes() == G_out.number_of_nodes()
        assert G_in.number_of_edges() == G_out.number_of_edges()
