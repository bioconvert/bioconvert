from bioconvert.fastq2fasta import Fastq2Fasta
from bioconvert import bioconvert_data
from bioconvert.core.compressor import make_in_gz_tester
from easydev import TempFile, md5
import pytest


# TODO: Add test of the unwrap_fasta method
@pytest.mark.parametrize("method", Fastq2Fasta.available_methods)
def test_conv(method):
    infile = bioconvert_data("test_fastq2fasta_v1.fastq")

    expected_outfile = bioconvert_data("test_fastq2fasta_v1.fasta")
    with TempFile(suffix=".fasta") as expected_unwrapped:
        Fastq2Fasta.unwrap_fasta(
            expected_outfile, expected_unwrapped.name, strip_comment=True)
        md5out = md5(expected_unwrapped.name)

    # One temporary file for the fasta created using the method
    # and one for an unwrapped version.
    # Some methods may output multi-line fasta, so we need to
    # compare md5 sums of unwrapped versions.
    with TempFile(suffix=".fasta") as outfile, \
            TempFile(suffix=".fasta") as unwrapped:
        convert = Fastq2Fasta(infile, outfile.name)
        convert(method=method)
        Fastq2Fasta.unwrap_fasta(
            outfile.name, unwrapped.name, strip_comment=True)
        assert md5(unwrapped.name) == md5out, \
            "{} failed".format(method)

@pytest.mark.parametrize(
    "method",
    filter(make_in_gz_tester(Fastq2Fasta), Fastq2Fasta.available_methods))
def test_in_gz(method):
    for sample_name in ["test_fastq2fasta_v1",
                        #"sample_v2", "sample_v3",  # GATB long headers bug
                        "sample_v4"]:
        infile = bioconvert_data("{}.fastq.gz".format(sample_name))

        expected_outfile = bioconvert_data("{}.fasta".format(sample_name))
        with TempFile(suffix=".fasta") as expected_unwrapped:
            Fastq2Fasta.unwrap_fasta(
                expected_outfile, expected_unwrapped.name, strip_comment=True)
            md5out = md5(expected_unwrapped.name)

        # One temporary file for the fasta created using the method
        # and one for an unwrapped version.
        # Some methods may output multi-line fasta, so we need to
        # compare md5 sums of unwrapped versions.
        with TempFile(suffix=".fasta") as outfile, \
                TempFile(suffix=".fasta") as unwrapped:
            convert = Fastq2Fasta(infile, outfile.name)
            convert(method=method)
            Fastq2Fasta.unwrap_fasta(
                outfile.name, unwrapped.name, strip_comment=True)
            assert md5(unwrapped.name) == md5out, \
                "{} failed for {}".format(method, sample_name)

@pytest.mark.parametrize("method", Fastq2Fasta.available_methods)
def _test_more_samples(method):
    for sample_name in ["sample_v2", "sample_v3", "sample_v4"]:
        infile = bioconvert_data("{}.fastq".format(sample_name))

        expected_outfile = bioconvert_data("{}.fasta".format(sample_name))
        with TempFile(suffix=".fasta") as expected_unwrapped:
            Fastq2Fasta.unwrap_fasta(
                expected_outfile, expected_unwrapped.name, strip_comment=True)
            md5out = md5(expected_unwrapped.name)

        # One temporary file for the fasta created using the method
        # and one for an unwrapped version.
        # Some methods may output multi-line fasta, so we need to
        # compare md5 sums of unwrapped versions.
        with TempFile(suffix=".fasta") as outfile, \
                TempFile(suffix=".fasta") as unwrapped:
            convert = Fastq2Fasta(infile, outfile.name)
            convert(method=method)
            Fastq2Fasta.unwrap_fasta(
                outfile.name, unwrapped.name, strip_comment=True)
            assert md5(unwrapped.name) == md5out, \
                "{} failed for {}".format(method, sample_name)
