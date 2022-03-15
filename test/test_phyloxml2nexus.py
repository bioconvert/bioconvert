import pytest
from easydev import TempFile, md5

from bioconvert.phyloxml2nexus import PHYLOXML2NEXUS

from . import test_dir

@pytest.mark.parametrize("method", PHYLOXML2NEXUS.available_methods)
def test_xml2nx_biopython(method):
    infile = f"{test_dir}/data/xml/{method}.xml"
    outfile = f"{test_dir}/data/nexus/{method}.nexus"
    with TempFile(suffix=".nexus") as tempfile:
        converter = PHYLOXML2NEXUS(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)