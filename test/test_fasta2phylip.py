import pytest
from easydev import TempFile, md5

from bioconvert.fasta2phylip import FASTA2PHYLIP

from . import test_dir


@pytest.mark.parametrize("method", FASTA2PHYLIP.available_methods)
def test_fa2phy_biopython(method):
    infile = f"{test_dir}/data/fasta/{method}.fasta"
    outfile = f"{test_dir}/data/bcf/{method}.phylip"
    with TempFile(suffix=".phylip") as tempfile:
        converter = FASTA2PHYLIP(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        # assert md5(tempfile.name) == md5(outfile)
        # md5 is difficult here https://github.com/bioconvert/bioconvert/issues/149
