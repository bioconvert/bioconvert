import os
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.phylip2stockholm import STOCKHOLM2PHYLIP, PHYLIP2STOCKHOLM

skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ
                                and os.environ['TRAVIS_PYTHON_VERSION'].startswith("2"), reason="On travis")


@skiptravis
def test_stockholm2phylip_biopython():
    infile = bioconvert_data("biopython.stockholm")
    outfile = bioconvert_data("biopython.phylip")
    with TempFile(suffix=".phylip") as tempfile:
        converter = STOCKHOLM2PHYLIP(infile, tempfile.name)
        converter(method='biopython')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


@skiptravis
def test_phylip2stockholm_biopython():
    infile = bioconvert_data("biopython.phylip")
    outfile = bioconvert_data("biopython.stockholm")
    with TempFile(suffix=".fasta") as tempfile:
        converter = PHYLIP2STOCKHOLM(infile, tempfile.name)
        converter(method='biopython')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)

@skiptravis
def test_stockholm2phylip_squizz():
    infile = bioconvert_data("squizz.stockholm")
    outfile = bioconvert_data("squizz.phylip")
    with TempFile(suffix=".phylip") as tempfile:
        converter = STOCKHOLM2PHYLIP(infile, tempfile.name)
        converter(method='squizz')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


@skiptravis
def test_phylip2stockholm_squizz():
    infile = bioconvert_data("squizz.phylip")
    outfile = bioconvert_data("squizz.stockholm")
    with TempFile(suffix=".stockholm") as tempfile:
        converter = PHYLIP2STOCKHOLM(infile, tempfile.name)
        converter(method='squizz')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
