import os
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.nexus2newick import NEXUS2NEWICK
import pytest


@pytest.mark.parametrize("method", NEXUS2NEWICK.available_methods)
def test_nx2nw_biopython(method):
    infile = bioconvert_data(method+".nexus")
    outfile = bioconvert_data(method+".newick")
    md5sum = md5(outfile)
    with TempFile(suffix=".newick") as tempfile:
        converter = NEXUS2NEWICK(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        if method == "biopython":
            assert md5(tempfile.name) == "accb8235ec0c5e6124de3837d8f684a9"
        else:
            assert md5(tempfile.name) == md5sum


