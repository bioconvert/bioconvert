import os
import bioconvert
from bioconvert import cram2bam
from bioconvert.cram2bam import CRAM2BAM
from bioconvert import bioconvert_data
from easydev import TempFile, md5
import pytest
from mock import patch

from bioconvert import bioconvert_data
from bioconvert.cram2sam import CRAM2SAM

reference = bioconvert_data("test_measles.fa")

@patch('bioconvert.cram2sam.input', return_value=reference)
def test_conv(x):
    infile = bioconvert_data("test_measles.cram")
    outfile = bioconvert_data("test_measles.sam")

    with TempFile(suffix=".sam") as tempfile:
        convert = CRAM2SAM(infile, tempfile.name)
        convert(method="samtools", reference=reference)

    with TempFile(suffix=".sam") as tempfile:
        convert = CRAM2SAM(infile, tempfile.name)
        convert(method="samtools")

@patch('bioconvert.cram2sam.input', return_value="not_found")
def test_conv_error(x):
    infile = bioconvert_data("test_measles.cram")
    outfile = bioconvert_data("test_measles.sam")
    with TempFile(suffix=".sam") as tempfile:
        convert = CRAM2SAM(infile, tempfile.name)
        try:
            convert(method="samtools")
            assert 0
        except IOError:
            assert 1


