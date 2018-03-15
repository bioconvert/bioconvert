import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.fasta2phylip import FASTA2PHYLIP


@pytest.mark.parametrize("method", FASTA2PHYLIP.available_methods)
def test_fa2phy_biopython(method):
    infile = bioconvert_data(method + ".fasta")
    outfile = bioconvert_data(method + ".phylip")
    with TempFile(suffix=".phylip") as tempfile:
        converter = FASTA2PHYLIP(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
