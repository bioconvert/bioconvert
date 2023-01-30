import pytest

from bioconvert.bam2fastq import BAM2FASTQ
from bioconvert import TempFile, md5

from . import test_dir


@pytest.mark.parametrize("method", BAM2FASTQ.available_methods)
def test_conv_all_methods(method):

    infile = f"{test_dir}/data/bam/test_measles.sorted.bam"
    with TempFile(suffix=".fastq") as tempfile:
        convert = BAM2FASTQ(infile, tempfile.name)
        convert(method=method)

    for ext in ["gz", "bz2", "dsrc"]:
        with TempFile(suffix=".fastq.{}".format(ext)) as tempfile:
            convert = BAM2FASTQ(infile, tempfile.name)
            convert(method=method)

    infile = f"{test_dir}/data/bam/test_measles_unpaired.sorted.bam"
    with TempFile(suffix=".fastq") as tempfile:
        convert = BAM2FASTQ(infile, tempfile.name)
        convert(method=method)

    for ext in ["gz", "bz2", "dsrc"]:
        with TempFile(suffix=".fastq.{}".format(ext)) as tempfile:
            convert = BAM2FASTQ(infile, tempfile.name)
            convert(method=method)


def test_method_bedtools():

    infile = f"{test_dir}/data/bam/test_measles.sorted.bam"
    with TempFile(suffix=".fastq") as tempfile:
        convert = BAM2FASTQ(infile, tempfile.name)
        convert(method="bedtools")


test_method_bedtools()
