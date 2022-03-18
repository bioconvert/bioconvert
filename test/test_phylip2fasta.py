import pytest
from easydev import TempFile, md5

from bioconvert.phylip2fasta import PHYLIP2FASTA

from . import test_dir

@pytest.mark.parametrize("method", PHYLIP2FASTA.available_methods)
def test_phy2pfa_biopython(method):
    infile = f"{test_dir}/data/phylip/biopython.phylip"
    outfile = f"{test_dir}/data/fasta/biopython.fasta"
    with TempFile(suffix=".phylip") as tempfile:
        converter = PHYLIP2FASTA(infile, tempfile.name)
        converter(method=method)
         # Check that the output is correct with a checksum
        if method in ["squizz", "goalign"]:
            assert md5(tempfile.name) == "98b661ab0f22ab6d4172135686cfc363"
        else:
            assert md5(tempfile.name) == md5(outfile)