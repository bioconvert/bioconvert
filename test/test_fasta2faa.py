from bioconvert.fasta2faa import FASTA2FAA
from bioconvert import bioconvert_data
from easydev import TempFile, md5
import pytest


def test_conv():
    infile = bioconvert_data("test_fasta2faa.fasta")
    expected_outfile = bioconvert_data("test_fasta2faa.faa")

    with TempFile(suffix=".faa") as outfile:
        convert = FASTA2FAA(infile, outfile.name)
        convert(method="bioconvert")
        assert md5(outfile.name) == md5(expected_outfile)

