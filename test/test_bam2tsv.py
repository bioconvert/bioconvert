import pytest

from bioconvert.bam2tsv import BAM2TSV
from bioconvert import bioconvert_data
from easydev import TempFile, md5


@pytest.mark.skipif(len(BAM2TSV.available_methods) == 0, reason="missing dependencies")
def test_bam2tsv():
    #your code here
    # you will need data for instance "mydata.fastq and mydata.fasta".
    # Put it in bioconvert/bioconvert/data
    # you can then use ::
    infile = bioconvert_data("test_measles.sorted.bam")
    #expected_outfile = bioconvert_data("test_measles.tsv")
    with TempFile(suffix=".tsv") as tempfile:
        convert = BAM2TSV(infile, tempfile.name)
        # Check that the output is correct with a checksum
        #assert md5(tempfile.name) == md5(expected_outfile)
        convert(method="pysam")
        assert md5(tempfile.name) == "4c5f3336be8a03c95a6c56be28581fb7"
        convert(method="samtools")
        assert md5(tempfile.name) == "4c5f3336be8a03c95a6c56be28581fb7"

