import pytest

from bioconvert.bam2fastq import BAM2FASTQ
from bioconvert import bioconvert_data
from easydev import TempFile, md5


@pytest.mark.skipif(len(BAM2FASTQ.available_methods) == 0, reason="missing dependencies")
def test_conv():
    infile = bioconvert_data("test_measles.sorted.bam")
    with TempFile(suffix=".fq") as tempfile:
        convert = BAM2FASTQ(infile, tempfile.name)
        convert()

        # Check that the output is correct with a checksum
        # Note that we cannot test the md5 on a gzip file but only 
        # on the original data. This check sum was computed
        # fro the unzipped version of biokit/data/converters/measles.bed
        assert md5(tempfile.name) == "8683ad696e52e4af67670d2631af6d1f"


@pytest.mark.parametrize("method", BAM2FASTQ.available_methods)
def test_conv_all_methods(method):
    infile = bioconvert_data("test_measles.sorted.bam")
    with TempFile(suffix=".cram") as tempfile:
        convert = BAM2FASTQ(infile, tempfile.name)
        convert(method=method)
