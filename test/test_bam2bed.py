from biokit.converters.bam2bed import Bam2Bed
from biokit import biokit_data
from easydev import TempFile, md5

def test_conv():
    infile = biokit_data("converters/measles.sorted.bam")
    with TempFile(suffix=".bed") as tempfile:
        convert = Bam2Bed(infile, tempfile.name)
        convert()

        # Check that the output is correct with a checksum
        # Note that we cannot test the md5 on a gzip file but only 
        # on the original data. This check sum was computed
        # fro the unzipped version of biokit/data/converters/measles.bed
        assert md5(tempfile.name) == "84702e19ba3a27900f271990e0cc72a0"
