import os
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.phylip2nexus import PHYLIP2NEXUS
import pytest


@pytest.mark.parametrize("method", PHYLIP2NEXUS.available_methods)
def test_phy2nx_biopython(method):
    infile = bioconvert_data(method+".phylip")
    outfile = bioconvert_data(method+".nexus")
    with TempFile(suffix=".nexus") as tempfile:
        converter = PHYLIP2NEXUS(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
