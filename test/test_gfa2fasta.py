from bioconvert.gfa2fasta import GFA2FASTA
from bioconvert import bioconvert_data
from easydev import TempFile, md5
import pytest

@pytest.mark.parametrize("method", GFA2FASTA.available_methods)
def test_conv_1(method):
    infile = bioconvert_data("test_gfa2fasta_v1.gfa")
    with TempFile(suffix=".fasta") as tempfile:
        convert = GFA2FASTA(infile, tempfile.name)
        convert()

        # Check that the output is correct with a checksum
        # Note that we cannot test the md5 on a gzip file but only 
        # on the original data. 
        assert md5(tempfile.name) == "6a02cf3308bf8d177a7d4cc55f317782"



@pytest.mark.parametrize("method", GFA2FASTA.available_methods)
def test_conv(method):
    infile = bioconvert_data("test_gfa2fasta.gfa")
    expected_outfile = bioconvert_data("test_gfa2fasta.fasta")
    md5out = md5(expected_outfile)

    # One temporary file for the fasta created using the method
    # and one for an unwrapped version.
    # Some methods may output multi-line fasta, so we need to
    # compare md5 sums of unwrapped versions.
    with TempFile(suffix=".fasta") as outfile:
        convert = GFA2FASTA(infile, outfile.name)
        convert(method=method)
        assert md5(convert.outfile) == md5out, \
            "{} failed".format(method)
