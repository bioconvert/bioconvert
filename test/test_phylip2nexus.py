import pytest
from bioconvert import TempFile, md5

from bioconvert.phylip2nexus import PHYLIP2NEXUS

from . import test_dir


@pytest.mark.parametrize("method", PHYLIP2NEXUS.available_methods)
def test_phy2nx_biopython(method):
    infile = f"{test_dir}/data/phylip/{method}.phylip"
    outfile = f"{test_dir}/data/nexus/{method}.nexus"
    with TempFile(suffix=".nexus") as tempfile:
        converter = PHYLIP2NEXUS(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
