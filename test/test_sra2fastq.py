import os
import pytest
import gzip
import shutil

from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.sra2fastq import SRA2FASTQ

@pytest.mark.parametrize("method", SRA2FASTQ.available_methods)
def test_sra2fastq_gz(method):
    infile = "SRR390728"
    outfile = bioconvert_data("SRR390728_1.fastq")
    outfile2 = bioconvert_data("SRR390728_2.fastq")
    with TempFile(suffix=".fastq.gz") as tempfile:
        converter = SRA2FASTQ(infile, tempfile.name, True)
        converter(method=method)
        outbasename,ext=os.path.splitext(tempfile.name)
        if ext == ".gz":
            outbasename,ext=os.path.splitext(outbasename)

        with gzip.open(outbasename+"_1.fastq.gz", 'rb') as f_in, open(outbasename+"_1.fastq", 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        with gzip.open(outbasename+"_2.fastq.gz", 'rb') as f_in, open(outbasename+"_2.fastq", 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
            
        # Check that the output is correct with a checksum
        assert md5(outbasename+"_1.fastq") == md5(outfile)
        assert md5(outbasename+"_2.fastq") == md5(outfile2)

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
    outfile = bioconvert_data("SRR6477205.fastq")
    
    with TempFile(suffix=".fastq.gz") as tempfile:
        converter = SRA2FASTQ(infile, tempfile.name, True)
        converter(method=method)

        outbasename=os.path.splitext(tempfile.name)[0]
        with gzip.open(tempfile.name, 'rb') as f_in, open(outbasename+".fastq", 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        
        # Check that the output is correct with a checksum
        assert md5(outbasename+".fastq") == md5(outfile)

@pytest.mark.parametrize("method", SRA2FASTQ.available_methods)
def test_sra2fastq_single(method):
    infile = "SRR6477205"
    outfile = bioconvert_data("SRR6477205.fastq")
    with TempFile(suffix=".fastq") as tempfile:
        converter = SRA2FASTQ(infile, tempfile.name, True)
        converter(method=method)
        
        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


