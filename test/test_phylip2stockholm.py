import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.phylip2stockholm import PHYLIP2STOCKHOLM


def test_phylip2stockholm_biopython():
    infile = bioconvert_data("biopython.phylip")
    outfile = bioconvert_data("biopython.stockholm")
    with TempFile(suffix=".fasta") as tempfile:
        converter = PHYLIP2STOCKHOLM(infile, tempfile.name)
        converter(method='biopython')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


@pytest.mark.skipif(PHYLIP2STOCKHOLM._method_squizz.is_disabled, reason="missing dependencies")
def test_phylip2stockholm_squizz():
    infile = bioconvert_data("squizz.phylip")
    outfile = bioconvert_data("squizz.stockholm")
    with TempFile(suffix=".stockholm") as tempfile:
        converter = PHYLIP2STOCKHOLM(infile, tempfile.name)
        converter(method='squizz')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
