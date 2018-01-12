import os
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.nexus2phylip import NEXUS2PHYLIP
import pytest


#skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ
#   and os.environ['TRAVIS_PYTHON_VERSION'].startswith("2"), reason="On travis")


#@skiptravis
@pytest.mark.parametrize("method", NEXUS2PHYLIP.available_methods)
def test_nx2phy_biopython(method):
    infile = bioconvert_data("goalign.nx")
    outfile = bioconvert_data("goalign.phy")
    with TempFile(suffix=".phy") as tempfile:
        converter = NEXUS2PHYLIP(infile, tempfile.name)
        converter(method='goalign')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


