import os
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.stockholm2phylip import STOCKHOLM2PHYLIP

#skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ
#                                and os.environ['TRAVIS_PYTHON_VERSION'].startswith("2"), reason="On travis")
#

#@skiptravis
def test_stockholm2phylip_biopython():
    infile = bioconvert_data("biopython.stockholm")
    outfile = bioconvert_data("biopython.phylip")
    with TempFile(suffix=".phylip") as tempfile:
        converter = STOCKHOLM2PHYLIP(infile, tempfile.name)
        converter(method='biopython')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


#@skiptravis
def test_stockholm2phylip_squizz():
    infile = bioconvert_data("squizz.stockholm")
    outfile = bioconvert_data("squizz.phylip")
    with TempFile(suffix=".phylip") as tempfile:
        converter = STOCKHOLM2PHYLIP(infile, tempfile.name)
        converter(method='squizz')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)

