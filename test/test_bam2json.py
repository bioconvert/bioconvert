from bioconvert.bam2json import BAM2JSON
from bioconvert import bioconvert_data
from easydev import TempFile, md5

def test_conv():
    infile = bioconvert_data("test_measles.sorted.bam")
    outfile = bioconvert_data("test_measles.sorted.json")
    with TempFile(suffix=".json") as tempfile:
        convert = BAM2JSON(infile, tempfile.name)
        convert(method="bamtools")

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
