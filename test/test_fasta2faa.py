from bioconvert.fasta2faa import FASTA2FAA
from easydev import TempFile, md5
import pytest

from . import test_dir

def test_conv():
    infile = f"{test_dir}/data/fasta/test_fasta2faa.fasta"
    expected_outfile = f"{test_dir}/data/faa/test_fasta2faa.faa"

    with TempFile(suffix=".faa") as outfile:
        convert = FASTA2FAA(infile, outfile.name)
        convert(method="bioconvert")
        assert md5(outfile.name) == md5(expected_outfile)

