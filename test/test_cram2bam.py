import os
from bioconvert.cram2bam import CRAM2BAM
from bioconvert import bioconvert_data
from easydev import TempFile, md5
import pytest



@pytest.mark.parametrize("method", CRAM2BAM.available_methods)
def test_conv(method):
    infile = bioconvert_data("test_measles.cram")
    outfile = bioconvert_data("test_measles.bam")
    #reference = bioconvert_data("test_measles.fa")

    with TempFile(suffix=".bam") as tempfile:
        convert = CRAM2BAM(infile, tempfile.name)
        convert(method=method)
