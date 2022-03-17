import pytest

from bioconvert.bam2sam import BAM2SAM
from easydev import TempFile, md5

from . import test_dir

# FIXME fails on CI action (no sequence)
@pytest.mark.skipif(len(BAM2SAM.available_methods) == 0, reason="missing dependencies")
def test_conv():
    infile = f"{test_dir}/data/bam/test_measles.sorted.bam"
    #outfile = biokit_data("converters/measles.sam")
    with TempFile(suffix=".bam") as tempfile:
        convert = BAM2SAM(infile, tempfile.name)
        # Test samtools
        convert(method="samtools")

        # Check that the output is correct with a checksum
        # Note that we cannot test the md5 on a gzip file but only 
        # on the original data. This check sum was computed
        # fro the unzipped version of biokit/data/converters/measles.bed
        #assert md5(tempfile.name) == md5(outfile)
        # output is a SAM that can be read and must have 
        import pysam
        sam = pysam.AlignmentFile(tempfile.name)
        assert sam.count() == 60

        # Test pysam 
        convert(method="pysam")

        convert = BAM2SAM(infile, tempfile.name)

        # Test sambamba
        convert(method="sambamba")
        # assert md5(tempfile.name) == "ad83af4d159005a77914c5503bc43802"