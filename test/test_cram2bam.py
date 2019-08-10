import os
import bioconvert
from bioconvert import cram2bam
from bioconvert.cram2bam import CRAM2BAM
from bioconvert import bioconvert_data
from easydev import TempFile, md5
import pytest
from mock import patch


reference = bioconvert_data("test_measles.fa")

@patch('bioconvert.cram2bam.input', return_value=reference)
def test_conv(x):
    infile = bioconvert_data("test_measles.bam")
    outfile = bioconvert_data("test_measles.cram")

    with TempFile(suffix=".bam") as tempfile:
        convert = CRAM2BAM(infile, tempfile.name)
        convert(method="samtools", reference=reference)

    with TempFile(suffix=".bam") as tempfile:
        convert = CRAM2BAM(infile, tempfile.name)
        convert(method="samtools")

@patch('bioconvert.cram2bam.input', return_value="not_found")
def test_conv_error(x):
    infile = bioconvert_data("test_measles.cram")
    outfile = bioconvert_data("test_measles.bam")
    with TempFile(suffix=".bam") as tempfile:
        convert = CRAM2BAM(infile, tempfile.name)
        try:
            convert(method="samtools")
            assert 0
        except IOError:
            assert 1


