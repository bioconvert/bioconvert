import os
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.phylip2clustal import PHYLIP2CLUSTAL

def test_phylip2clustal_biopython():
    infile = bioconvert_data("biopython.phylip")
    outfile = bioconvert_data("biopython.clustal")
    with TempFile(suffix=".fasta") as tempfile:
        converter = PHYLIP2CLUSTAL(infile, tempfile.name)
        converter(method='biopython')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


@pytest.mark.skipif(PHYLIP2CLUSTAL._method_squizz.is_disabled, reason="missing dependencies")
def test_phylip2clustal_squizz():
    infile = bioconvert_data("squizz.phylip")
    outfile = bioconvert_data("squizz.clustal")
    with TempFile(suffix=".clustal") as tempfile:
        converter = PHYLIP2CLUSTAL(infile, tempfile.name)
        converter(method='squizz')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
