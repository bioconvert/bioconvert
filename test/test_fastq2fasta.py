from bioconvert.fastq2fasta import FASTQ2FASTA
from bioconvert.core.decorators import make_in_gz_tester, requires
from easydev import TempFile, md5
import pytest

from . import test_dir

# TODO: Add test of the unwrap_fasta method
@pytest.mark.parametrize("method", FASTQ2FASTA.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/fastq2fasta/test_fastq2fasta_v1.fastq"

    expected_outfile = f"{test_dir}/data/fastq2fasta/test_fastq2fasta_v1.fasta"
    with TempFile(suffix=".fasta") as expected_unwrapped:
        FASTQ2FASTA.unwrap_fasta(
            expected_outfile, expected_unwrapped.name, strip_comment=True)
        md5out = md5(expected_unwrapped.name)

    # One temporary file for the fasta created using the method
    # and one for an unwrapped version.
    # Some methods may output multi-line fasta, so we need to
    # compare md5 sums of unwrapped versions.
    with TempFile(suffix=".fasta") as outfile, \
            TempFile(suffix=".fasta") as unwrapped:
        convert = FASTQ2FASTA(infile, outfile.name)
        convert(method=method)
        FASTQ2FASTA.unwrap_fasta(
            outfile.name, unwrapped.name, strip_comment=True)
        assert md5(unwrapped.name) == md5out, \
            "{} failed".format(method)




@pytest.mark.parametrize(
    "method",
    filter(make_in_gz_tester(FASTQ2FASTA), FASTQ2FASTA.available_methods))
@pytest.mark.skipif(
    requires(external_binary="unpigz")(object()).is_disabled,
    reason="missing dependencies",
)
def test_in_gz(method):
    for sample_name in ["test_fastq2fasta_v1",
                        "sample_v2", "sample_v3",  
                        "sample_v4"]:
        infile = f"{test_dir}/data/fastq2fasta/{sample_name}.fastq.gz"

        expected_outfile = f"{test_dir}/data/fastq2fasta/{sample_name}.fasta"
        with TempFile(suffix=".fasta") as expected_unwrapped:
            FASTQ2FASTA.unwrap_fasta(
                expected_outfile, expected_unwrapped.name, strip_comment=True)
            md5out = md5(expected_unwrapped.name)

        # One temporary file for the fasta created using the method
        # and one for an unwrapped version.
        # Some methods may output multi-line fasta, so we need to
        # compare md5 sums of unwrapped versions.
        with TempFile(suffix=".fasta") as outfile, \
                TempFile(suffix=".fasta") as unwrapped:
            convert = FASTQ2FASTA(infile, outfile.name)
            convert(method=method)
            FASTQ2FASTA.unwrap_fasta(
                outfile.name, unwrapped.name, strip_comment=True)
            assert md5(unwrapped.name) == md5out, \
                "{} failed for {}".format(method, sample_name)

@pytest.mark.parametrize("method", FASTQ2FASTA.available_methods)
def test_more_samples(method):
    for sample_name in ["sample_v2", "sample_v3", "sample_v4"]:

        infile = f"{test_dir}/data/fastq2fasta/{sample_name}.fastq"

        expected_outfile = f"{test_dir}/data/fastq2fasta/{sample_name}.fasta"
        with TempFile(suffix=".fasta") as expected_unwrapped:
            FASTQ2FASTA.unwrap_fasta(
                expected_outfile, expected_unwrapped.name, strip_comment=True)
            md5out = md5(expected_unwrapped.name)

        # One temporary file for the fasta created using the method
        # and one for an unwrapped version.
        # Some methods may output multi-line fasta, so we need to
        # compare md5 sums of unwrapped versions.
        with TempFile(suffix=".fasta") as outfile, \
                TempFile(suffix=".fasta") as unwrapped:
            convert = FASTQ2FASTA(infile, outfile.name)
            convert(method=method)
            FASTQ2FASTA.unwrap_fasta(
                outfile.name, unwrapped.name, strip_comment=True)
            assert md5(unwrapped.name) == md5out, \
                "{} failed for {}".format(method, sample_name)
