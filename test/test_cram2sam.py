import os
from bioconvert.cram2sam import CRAM2SAM
from bioconvert import bioconvert_data
from easydev import TempFile, md5


def test_conv():
    infile = bioconvert_data("test_measles.cram")
    outfile = bioconvert_data("test_measles.sam")
    reference = bioconvert_data("test_measles.fa")

    with TempFile(suffix=".sam") as tempfile:
        convert = CRAM2SAM(infile, tempfile.name, reference)
        convert(method="samtools")
