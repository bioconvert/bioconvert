from bioconvert.gfa2fasta import GFA2FASTA
from easydev import TempFile, md5
import pytest

from . import test_dir

@pytest.mark.parametrize("method", GFA2FASTA.available_methods)
def test_conv_1(method):
    infile = f"{test_dir}/data/gfa/test_gfa2fasta_v1.gfa"
    with TempFile(suffix=".fasta") as tempfile:
        convert = GFA2FASTA(infile, tempfile.name)
        convert()

        # Check that the output is correct with a checksum
        # Note that awk wrap fasta while python method does not 
        assert md5(tempfile.name) in \
            [ "6a02cf3308bf8d177a7d4cc55f317782",
              'aced24d965c5420781ad0884a56ed3c5']




@pytest.mark.parametrize("method", GFA2FASTA.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/gfa/test_gfa2fasta.gfa"
    expected_outfile = f"{test_dir}/data/fasta/test_gfa2fasta.fasta"
    md5out = md5(expected_outfile)

    with TempFile(suffix=".fasta") as outfile:
        convert = GFA2FASTA(infile, outfile.name)
        convert(method=method)
        assert md5(convert.outfile) == md5out, \
            "{} failed".format(method)