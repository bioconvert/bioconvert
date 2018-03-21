import os
from bioconvert.bam2cram import BAM2CRAM
from bioconvert import bioconvert_data
from easydev import TempFile, md5
import pytest



@pytest.mark.parametrize("method", BAM2CRAM.available_methods)
def test_conv(method):
    infile = bioconvert_data("test_measles.bam")
    outfile = bioconvert_data("test_measles.cram")
    #reference = bioconvert_data("test_measles.fa")

    with TempFile(suffix=".cram") as tempfile:
        convert = BAM2CRAM(infile, tempfile.name)
        convert(method=method)
