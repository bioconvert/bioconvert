import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.newick2nexus import NEWICK2NEXUS


@pytest.mark.parametrize("method", NEWICK2NEXUS.available_methods)
def test_nw2nx_biopython(method):
    infile = bioconvert_data(method + ".newick")
    outfile = bioconvert_data(method + ".nexus")
    with TempFile(suffix=".nexus") as tempfile:
        converter = NEWICK2NEXUS(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
