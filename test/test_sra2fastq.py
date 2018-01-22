import os
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.sra2fastq import SRA2FASTQ
import pytest

@pytest.mark.parametrize("method", SRA2FASTQ.available_methods)
def test_sra2fastq_gz(method):
    infile = "SRR390728"
    outfile = bioconvert_data("SRR390728_1.fastq.gz")
    outfile2 = bioconvert_data("SRR390728_2.fastq.gz")
    with TempFile(suffix=".fastq.gz") as tempfile:
        converter = SRA2FASTQ(infile, tempfile.name, True)
        converter(method=method)
        outbasename,ext=os.path.splitext(tempfile.name)
        if ext == ".gz":
            outbasename,ext=os.path.splitext(outbasename)
        
        # Check that the output is correct with a checksum
        assert md5(outbasename+"_1.fastq.gz") == md5(outfile)
        assert md5(outbasename+"_2.fastq.gz") == md5(outfile2)

@pytest.mark.parametrize("method", SRA2FASTQ.available_methods)
def test_sra2fastq(method):
    infile = "SRR390728"
    outfile = bioconvert_data("SRR390728_1.fastq")
    outfile2 = bioconvert_data("SRR390728_2.fastq")
    with TempFile(suffix=".fastq") as tempfile:
        converter = SRA2FASTQ(infile, tempfile.name, True)
        converter(method=method)
        outbasename=os.path.splitext(tempfile.name)[0]
        
        # Check that the output is correct with a checksum
        assert md5(outbasename+"_1.fastq") == md5(outfile)
        assert md5(outbasename+"_2.fastq") == md5(outfile2)

@pytest.mark.parametrize("method", SRA2FASTQ.available_methods)
def test_sra2fastq_gz_single(method):
    infile = "SRR6477205"
    outfile = bioconvert_data("SRR6477205.fastq.gz")
    with TempFile(suffix=".fastq.gz") as tempfile:
        converter = SRA2FASTQ(infile, tempfile.name, True)
        converter(method=method)
        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)

@pytest.mark.parametrize("method", SRA2FASTQ.available_methods)
def test_sra2fastq_single(method):
    infile = "SRR6477205"
    outfile = bioconvert_data("SRR6477205.fastq")
    with TempFile(suffix=".fastq") as tempfile:
        converter = SRA2FASTQ(infile, tempfile.name, True)
        converter(method=method)
        
        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


