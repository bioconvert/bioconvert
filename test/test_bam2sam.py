from bioconvert.bam2sam import BAM2SAM
from bioconvert import bioconvert_data
from easydev import TempFile, md5


def test_conv():
    infile = bioconvert_data("test_measles.sorted.bam")
    #outfile = biokit_data("converters/measles.sam")
    with TempFile(suffix=".bam") as tempfile:
        convert = BAM2SAM(infile, tempfile.name)
        convert()

        # Check that the output is correct with a checksum
        # Note that we cannot test the md5 on a gzip file but only 
        # on the original data. This check sum was computed
        # fro the unzipped version of biokit/data/converters/measles.bed
        #assert md5(tempfile.name) == md5(outfile)
        # output is a SAM that can be read and must have 
        import pysam
        sam = pysam.AlignmentFile(tempfile.name)
        assert sam.count() == 60
