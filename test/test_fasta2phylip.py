import os
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.fasta2phylip import FASTA2PHYLIP, PHYLIP2FASTA

skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ
                                and os.environ['TRAVIS_PYTHON_VERSION'].startswith("2"), reason="On travis")


@skiptravis
def test_fa2phy_biopython():
    infile = bioconvert_data("biopython.fasta")
    outfile = bioconvert_data("biopython.phylip")
    with TempFile(suffix=".phylip") as tempfile:
        converter = FASTA2PHYLIP(infile, tempfile.name)
        converter(method='biopython')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


@skiptravis
def test_phy2fa_biopython():
    infile = bioconvert_data("biopython.phylip")
    outfile = bioconvert_data("biopython.fasta")
    with TempFile(suffix=".fasta") as tempfile:
        converter = PHYLIP2FASTA(infile, tempfile.name)
        converter(method='biopython')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)

@skiptravis
def test_fa2phy_squizz():
    infile = bioconvert_data("squizz.fasta")
    outfile = bioconvert_data("squizz.phylip")
    with TempFile(suffix=".phylip") as tempfile:
        converter = FASTA2PHYLIP(infile, tempfile.name)
        converter(method='squizz')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


@skiptravis
def test_phy2fa_squizz():
    infile = bioconvert_data("squizz.phylip")
    outfile = bioconvert_data("squizz.fasta")
    with TempFile(suffix=".fasta") as tempfile:
        converter = PHYLIP2FASTA(infile, tempfile.name)
        converter(method='squizz')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
