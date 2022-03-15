import os
import bioconvert
from bioconvert import bam2cram
from bioconvert.bam2cram import BAM2CRAM
from easydev import TempFile, md5
import pytest
from mock import patch

from . import test_dir

reference = f"{test_dir}/data/fasta/test_measles.fa"

# may not work with patch decorator. there is only one method anyway.
#@pytest.mark.parametrize("method", BAM2CRAM.available_methods)

@patch('bioconvert.bam2cram.input', return_value=reference)
def test_conv(x):
    infile = f"{test_dir}/data/bam/test_measles.bam"
    outfile = f"{test_dir}/data/cram/test_measles.cram"

    with TempFile(suffix=".cram") as tempfile:
        convert = BAM2CRAM(infile, tempfile.name)
        convert(method="samtools", reference=reference)

    with TempFile(suffix=".cram") as tempfile:
        convert = BAM2CRAM(infile, tempfile.name)
        convert(method="samtools")

@patch('bioconvert.bam2cram.input', return_value="not_found")
def test_conv_error(x):
    infile = f"{test_dir}/data/bam/test_measles.bam"
    outfile = f"{test_dir}/data/cram/test_measles.cram"
    with TempFile(suffix=".cram") as tempfile:
        convert = BAM2CRAM(infile, tempfile.name)
        try:
            convert(method="samtools")
            assert 0
        except IOError:
            assert 1


