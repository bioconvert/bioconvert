from bioconvert.bam2fastq import BAM2Fastq
from bioconvert import bioconvert_data
from easydev import TempFile, md5
import pytest

@pytest.mark.parametrize("method", BAM2Fastq.available_methods)
def test_conv(method):
    infile = bioconvert_data("test_measles.sorted.bam")
    with TempFile(suffix=".fq") as tempfile:
        convert = BAM2Fastq(infile, tempfile.name)
        convert(method=method)

        # Check that the output is correct with a checksum
        # Note that we cannot test the md5 on a gzip file but only 
        # on the original data. This check sum was computed
        # fro the unzipped version of biokit/data/converters/measles.bed
        assert md5(tempfile.name) == "8683ad696e52e4af67670d2631af6d1f"


        # for method in convert.available_methods:
        #     convert(method=method)
