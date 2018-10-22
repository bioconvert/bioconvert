import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.sam2bam import SAM2BAM


@pytest.mark.skipif(len(SAM2BAM.available_methods) == 0, reason="missing dependencies")
def test_conv():
    infile = bioconvert_data("test_measles.sam")
    outfile = bioconvert_data("test_measles.bam")
    with TempFile(suffix=".bam") as tempfile:
        convert = SAM2BAM(infile, tempfile.name)
        convert()

        # Check that the output is correct with a checksum
        # Note that we cannot test the md5 on a gzip file but only 
        # on the original data. This check sum was computed
        # fro the unzipped version of biokit/data/converters/measles.bed
        #assert md5(tempfile.name) == md5(outfile)
        # 5cd453e698bccf942431618c945c226e

        # TODO FIXME this md5sum changes with time probably due to the samtools
        # version being encoded in the header. Need to find another to 
        # check the content of the BAM/SAM file. For isntance using
        # bamtools/samtools stats or pysam. 

