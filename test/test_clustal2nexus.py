import os
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.clustal2nexus import CLUSTAL2NEXUS

skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ
                                and os.environ['TRAVIS_PYTHON_VERSION'].startswith("2"), reason="On travis")



@skiptravis
def test_clustal2nexus_goalign():
    infile = bioconvert_data("goalign.clustal")
    outfile = bioconvert_data("goalign.nexus")
    with TempFile(suffix=".nexus") as tempfile:
        converter = CLUSTAL2NEXUS(infile, tempfile.name)
        converter(method='goalign')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


