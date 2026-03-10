import os
import pytest
import gzip
import shutil

from bioconvert import TempFile, md5

from bioconvert.sra2fastq import SRA2FASTQ

from . import test_dir


@pytest.mark.slow
@pytest.mark.flaky
@pytest.mark.parametrize("method", SRA2FASTQ.available_methods)
def test_sra2fastq_gz(method):
    """
    :param method: the method to be tested among the available methods
    :type method: str

    .. note:: All output files are cut to 10 reads to reduce file size and speed up tests. So the test option = True
    """
    infile = "SRR37522531"
    with TempFile(suffix=".fastq.gz") as tempfile:
        # if test=True only the first 10 reads from the sra file will be taken
        converter = SRA2FASTQ(infile, tempfile.name, test=True)
        converter(method=method)
        outbasename, ext = os.path.splitext(tempfile.name)
        if ext == ".gz":
            outbasename, ext = os.path.splitext(outbasename)

        with gzip.open(outbasename + "_1.fastq.gz", "rb") as f_in, open(
            outbasename + "_1.fastq", "wb"
        ) as f_out:
            shutil.copyfileobj(f_in, f_out)
        with gzip.open(outbasename + "_2.fastq.gz", "rb") as f_in, open(
            outbasename + "_2.fastq", "wb"
        ) as f_out:
            shutil.copyfileobj(f_in, f_out)

        # Check that the output is correct with a checksum
        assert md5(outbasename + "_1.fastq") in ["8bb130ace420f7a4a11eec18c282c887","0a7fe8a6d4ac2a8c77fa00ef09d05305"]
        assert md5(outbasename + "_2.fastq") in ["dab86ad2c83ce7818613b61739bad960","81d42149953502a8efd5c7878dfdd1b9"]


@pytest.mark.slow
@pytest.mark.flaky
@pytest.mark.parametrize("method", SRA2FASTQ.available_methods)
def test_sra2fastq(method):
    infile = "SRR37522531"
    with TempFile(suffix=".fastq") as tempfile:
        converter = SRA2FASTQ(infile, tempfile.name, test=True)
        converter(method=method)
        outbasename = os.path.splitext(tempfile.name)[0]

        # Check that the output is correct with a checksum
        assert md5(outbasename + "_1.fastq") in ["8bb130ace420f7a4a11eec18c282c887", "0a7fe8a6d4ac2a8c77fa00ef09d05305"]
        assert md5(outbasename + "_2.fastq") in ["dab86ad2c83ce7818613b61739bad960", "81d42149953502a8efd5c7878dfdd1b9"]


@pytest.mark.slow
@pytest.mark.flaky
@pytest.mark.parametrize("method", SRA2FASTQ.available_methods)
def test_sra2fastq_gz_single(method):
    infile = "SRR30092023"

    with TempFile(suffix=".fastq.gz") as tempfile:
        converter = SRA2FASTQ(infile, tempfile.name, test=True)
        converter(method=method)

        outbasename = os.path.splitext(tempfile.name)[0]
        with gzip.open(tempfile.name, "rb") as f_in, open(
            outbasename + ".fastq", "wb"
        ) as f_out:
            shutil.copyfileobj(f_in, f_out)

        # Check that the output is correct with a checksum
        assert md5(outbasename + ".fastq") in [
            "8aa18792cfa90a43f4ffaa35cf290a5c",
            "06d0e45af223a0f25d96355a7e12012b"]

@pytest.mark.slow
@pytest.mark.flaky
@pytest.mark.parametrize("method", SRA2FASTQ.available_methods)
def test_sra2fastq_single(method):
    infile = "SRR30092023"
    with TempFile(suffix=".fastq") as tempfile:
        converter = SRA2FASTQ(infile, tempfile.name, test=True)
        converter(method=method)

        # Check that the output is correct with a checksum
        # fastq-dump uses only 10 reads hence a different md5sum
        assert md5(tempfile.name) in [
            "8aa18792cfa90a43f4ffaa35cf290a5c",
            "06d0e45af223a0f25d96355a7e12012b"]
