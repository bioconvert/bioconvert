from bioconvert.gfa2fasta import GFA2FASTA
from bioconvert import bioconvert_data
from easydev import TempFile, md5


def test_conv():
    infile = bioconvert_data("test_v1.gfa")
    with TempFile(suffix=".fasta") as tempfile:
        convert = GFA2FASTA(infile, tempfile.name)
        convert()

        # Check that the output is correct with a checksum
        # Note that we cannot test the md5 on a gzip file but only 
        # on the original data. 
        assert md5(tempfile.name) == "6a02cf3308bf8d177a7d4cc55f317782"
