from bioconvert.fasta2fastq import FASTA2FASTQ
from bioconvert import TempFile, md5
import pytest

from . import test_dir

# TODO: Add test of the unwrap_fasta method
@pytest.mark.parametrize("method", FASTA2FASTQ.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/fasta/test_fastq2fasta_v1.fasta"
    qual_file = f"{test_dir}/data/qual/test_fasta2fastq.qual"

    expected_outfile_no_qual = f"{test_dir}/data/fastq/test_fasta2fastq.fastq"
    md5out_no_qual = md5(expected_outfile_no_qual)

    expected_outfile_qual = f"{test_dir}/data/fastq/test_fastq2fasta_v1.fastq"
    md5out_qual = md5(expected_outfile_qual)

    # One temporary file for the fasta created using the method
    # and one for an unwrapped version.
    # Some methods may output multi-line fasta, so we need to
    # compare md5 sums of unwrapped versions.
    with TempFile(suffix=".fastq") as outfile:
        convert = FASTA2FASTQ(infile, outfile.name)
        convert(method=method)
        assert md5(outfile.name) == md5out_no_qual, "{} failed".format(method)

    with TempFile(suffix=".fastq") as outfile:
        convert = FASTA2FASTQ(infile, outfile.name)
        convert(method=method, quality_file=qual_file)
        assert md5(outfile.name) == md5out_qual, "{} failed".format(method)
