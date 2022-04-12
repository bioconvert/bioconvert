import os
import bioconvert
from bioconvert import cram2bam
from bioconvert.cram2bam import CRAM2BAM
from easydev import TempFile, md5
import pytest
from mock import patch

from . import test_dir

reference = f"{test_dir}/data/fasta/test_measles.fa"


@patch("bioconvert.cram2bam.input", return_value=reference)
def test_conv(x):
    infile = f"{test_dir}/data/bam/test_measles.bam"
    outfile = f"{test_dir}/data/cram/test_measles.cram"

    with TempFile(suffix=".bam") as tempfile:
        convert = CRAM2BAM(infile, tempfile.name)
        convert(method="samtools", reference=reference)

    with TempFile(suffix=".bam") as tempfile:
        convert = CRAM2BAM(infile, tempfile.name)
        convert(method="samtools")


@patch("bioconvert.cram2bam.input", return_value="not_found")
def test_conv_error(x):
    infile = f"{test_dir}/data/cram/test_measles.cram"
    outfile = f"{test_dir}/data/bam/test_measles.bam"
    with TempFile(suffix=".bam") as tempfile:
        convert = CRAM2BAM(infile, tempfile.name)
        try:
            convert(method="samtools")
            assert 0
        except IOError:
            assert 1
