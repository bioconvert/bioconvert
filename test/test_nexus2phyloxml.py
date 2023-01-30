import pytest
from bioconvert import TempFile, md5

from bioconvert.nexus2phyloxml import NEXUS2PHYLOXML

from . import test_dir


@pytest.mark.parametrize("method", NEXUS2PHYLOXML.available_methods)
def test_nx2xml_biopython(method):
    infile = f"{test_dir}/data/nexus/{method}.nexus"
    outfile = f"{test_dir}/data/phyloxml/{method}.xml"
    with TempFile(suffix=".phyloxml") as tempfile:
        converter = NEXUS2PHYLOXML(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
