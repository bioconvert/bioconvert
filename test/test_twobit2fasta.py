import os
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.twobit2fasta import TWOBIT2FASTA
import pytest


@pytest.mark.parametrize("method", TWOBIT2FASTA.available_methods)
def test_twobit2fasta_ucsc(method):
    infile = bioconvert_data("ucsc.2bit")
    outfile = bioconvert_data("ucsc.fasta")
    with TempFile(suffix=".fasta") as tempfile:
        converter = TWOBIT2FASTA(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
