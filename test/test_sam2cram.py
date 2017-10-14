import os
from bioconvert.sam2cram import SAM2CRAM
from bioconvert import bioconvert_data
from easydev import TempFile, md5


def test_conv():
    infile = bioconvert_data("test_measles.sam")
    outfile = bioconvert_data("test_measles.cram")
    reference = bioconvert_data("test_measles.fa")
    with TempFile(suffix=".cram") as tempfile:
        convert = SAM2CRAM(infile, tempfile.name, reference)
        convert()

        # Check that the output is correct with a checksum
        # Note that we cannot test the md5 on a gzip file but only 
        # on the original data. This check sum was computed
        # fro the unzipped version of bioconvert/data/converters/measles.bed
        #assert md5(tempfile.name) == md5(outfile)
        size = os.path.getsize(tempfile.name)
        # compressed file size may change. I have seen 6115, 6608, 6141
        assert size > 5800 and size < 7000

    with TempFile(suffix=".cram") as tempfile:
        convert = SAM2CRAM(infile, tempfile.name, reference=None)
        convert()

    with TempFile(suffix=".cram") as tempfile:
        convert = SAM2CRAM(infile, tempfile.name, reference="dummy.fa")
        convert()



