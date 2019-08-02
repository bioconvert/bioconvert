import pytest

from bioconvert.bam2fastq import BAM2FASTQ
from bioconvert import bioconvert_data
from easydev import TempFile, md5




@pytest.mark.parametrize("method", BAM2FASTQ.available_methods)
def test_conv_all_methods(method):
    infile = bioconvert_data("test_measles.sorted.bam")
    with TempFile(suffix=".fastq") as tempfile:
        convert = BAM2FASTQ(infile, tempfile.name)
        convert(method)

    for ext in ['gz', 'bz2', 'dsrc']:
        with TempFile(suffix=".fastq.{}".format(ext)) as tempfile:
            convert = BAM2FASTQ(infile, tempfile.name)
            convert(method)
