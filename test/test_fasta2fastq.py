from bioconvert.fasta2fastq import Fasta2Fastq
from bioconvert import bioconvert_data
from easydev import TempFile, md5
import pytest


# TODO: Add test of the unwrap_fasta method
@pytest.mark.parametrize("method", Fasta2Fastq.available_methods)
def test_conv(method):
    infile = bioconvert_data("test_fastq2fasta_v1.fasta")

    expected_outfile = bioconvert_data("test_fasta2fastq.fastq")
    md5out = md5(expected_outfile)

    # One temporary file for the fasta created using the method
    # and one for an unwrapped version.
    # Some methods may output multi-line fasta, so we need to
    # compare md5 sums of unwrapped versions.
    with TempFile(suffix=".fastq") as outfile:
        convert = Fasta2Fastq(infile, outfile.name)
        convert(method=method)
        assert md5(outfile.name) == md5out, \
            "{} failed".format(method)
