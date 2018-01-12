import os
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.phylip2fasta import PHYLIP2FASTA
import pytest


skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ
                                and os.environ['TRAVIS_PYTHON_VERSION'].startswith("2"), reason="On travis")


@skiptravis
@pytest.mark.parametrize("method", PHYLIP2FASTA.available_methods)
def test_fa2phy_biopython(method):
    infile = bioconvert_data("biopython.phylip")
    outfile = bioconvert_data("biopython.fasta")
    with TempFile(suffix=".phylip") as tempfile:
        converter = PHYLIP2FASTA(infile, tempfile.name)
        converter(method='biopython')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


