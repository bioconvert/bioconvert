import os
from bioconvert.sam2cram import SAM2CRAM
from easydev import TempFile, md5
import pytest
from mock import patch 

from . import test_dir

reference = f"{test_dir}/data/fasta/test_measles.fa"

@patch('bioconvert.sam2cram.input', return_value=reference)
def test_conv(x):
    infile = f"{test_dir}/data/cram/test_measles.cram"
    outfile = f"{test_dir}/data/sam/test_measles.sam"

    with TempFile(suffix=".cram") as tempfile:
        convert = SAM2CRAM(infile, tempfile.name, reference)
        convert(reference=reference)

        # Check that the output is correct with a checksum
        # Note that we cannot test the md5 on a gzip file but only
        # on the original data. This check sum was computed
        # fro the unzipped version of bioconvert/data/converters/measles.bed
        # assert md5(tempfile.name) == md5(outfile)
        size = os.path.getsize(tempfile.name)
        # compressed file size may change. I have seen 6115, 6608, 6141,9080
        assert size > 5800 and size < 10000


@patch('bioconvert.sam2cram.input', return_value="not_found")
def test_conv_error(x):
    infile = f"{test_dir}/data/cram/test_measles.cram"
    outfile = f"{test_dir}/data/sam/test_measles.sam"
    with TempFile(suffix=".sam") as tempfile:
        convert = SAM2CRAM(infile, tempfile.name)
        try:
            convert(method="samtools")
            assert 0
        except IOError:
            assert 1