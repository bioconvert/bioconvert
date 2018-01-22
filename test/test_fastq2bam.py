import os
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.fastq2bam import FASTQ2BAM
from bioconvert.core.compressor import make_in_gz_tester

@pytest.mark.parametrize("method", filter(make_in_gz_tester(FASTQ2BAM),FASTQ2BAM.available_methods))
def test_fastq2bam_paired_gz(method):
    infile  = bioconvert_data(method+"_1.fastq.gz")
    infile2 = bioconvert_data(method+"_2.fastq.gz")
    outfile = bioconvert_data(method+".bam")
    with TempFile(suffix=".bam") as tempfile:
        converter = FASTQ2BAM(infile, tempfile.name, infile2)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)

@pytest.mark.parametrize("method", filter(make_in_gz_tester(FASTQ2BAM),FASTQ2BAM.available_methods))
def test_fastq2bam_single_gz(method):
    infile  = bioconvert_data(method+"_1.fastq.gz")
    outfile = bioconvert_data(method+"_1.bam")
    with TempFile(suffix=".bam") as tempfile:
        converter = FASTQ2BAM(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)

@pytest.mark.parametrize("method", FASTQ2BAM.available_methods)
def test_fastq2bam_paired(method):
    infile  = bioconvert_data(method+"_1.fastq")
    infile2 = bioconvert_data(method+"_2.fastq")
    outfile = bioconvert_data(method+".bam")
    with TempFile(suffix=".bam") as tempfile:
        converter = FASTQ2BAM(infile, tempfile.name, infile2)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)

@pytest.mark.parametrize("method", FASTQ2BAM.available_methods)
def test_fastq2bam_single(method):
    infile  = bioconvert_data(method+"_1.fastq")
    outfile = bioconvert_data(method+"_1.bam")
    with TempFile(suffix=".bam") as tempfile:
        converter = FASTQ2BAM(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
