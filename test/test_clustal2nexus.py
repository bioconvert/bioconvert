import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.clustal2nexus import CLUSTAL2NEXUS


@pytest.mark.parametrize("method", CLUSTAL2NEXUS.available_methods)
def test_nx2aln(method):
    infile = bioconvert_data(method + ".clustal")
    outfile = bioconvert_data(method + ".nexus")
    with TempFile(suffix=".nexus") as tempfile:
        converter = CLUSTAL2NEXUS(infile, tempfile.name)
        converter(method=method)
        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
