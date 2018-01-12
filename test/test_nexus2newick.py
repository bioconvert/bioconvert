import os
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.nexus2newick import NEXUS2NEWICK
import pytest


skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ
                                and os.environ['TRAVIS_PYTHON_VERSION'].startswith("2"), reason="On travis")


@skiptravis
@pytest.mark.parametrize("method", NEXUS2NEWICK.available_methods)
def test_nx2nw_biopython(method):
    infile = bioconvert_data("gotree.nx")
    outfile = bioconvert_data("gotree.nw")
    with TempFile(suffix=".nw") as tempfile:
        converter = NEXUS2NEWICK(infile, tempfile.name)
        converter(method='gotree')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


