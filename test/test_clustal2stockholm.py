import os
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.clustal2stockholm import CLUSTAL2STOCKHOLM

skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ
                                and os.environ['TRAVIS_PYTHON_VERSION'].startswith("2"), reason="On travis")


@skiptravis
def test_clustal2stockholm_biopython():
    infile = bioconvert_data("biopython.clustal")
    outfile = bioconvert_data("biopython.stockholm")
    with TempFile(suffix=".stockholm") as tempfile:
        converter = CLUSTAL2STOCKHOLM(infile, tempfile.name)
        converter(method='biopython')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


@skiptravis
def test_clustal2stockholm_squizz():
    infile = bioconvert_data("squizz.clustal")
    outfile = bioconvert_data("squizz.stockholm")
    with TempFile(suffix=".stockholm") as tempfile:
        converter = CLUSTAL2STOCKHOLM(infile, tempfile.name)
        converter(method='squizz')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


