import os
from bioconvert.sam2paf import SAM2PAF
from bioconvert import bioconvert_data
from easydev import TempFile, md5

where = "testing/sam2paf"

def test_conv():
    infile = bioconvert_data("test_sam2paf_v1.sam", where)
    outfile = bioconvert_data("test_sam2paf_v1.paf", where)
    checksum = md5(outfile)

    with TempFile(suffix=".paf") as tempfile:
        convert = SAM2PAF(infile, tempfile.name)
        convert()
        assert checksum == md5(tempfile.name)
        assert convert.skipped == 17


def test_bad1_input():
    infile = bioconvert_data("test_sam2paf_bad1.sam", where)
    with TempFile(suffix=".paf") as tempfile:
        convert = SAM2PAF(infile, tempfile.name)
        try:
            convert()
            assert False
        except ValueError:
            assert True

def test_bad2_input():
    infile = bioconvert_data("test_sam2paf_bad2.sam", where)
    with TempFile(suffix=".paf") as tempfile:
        convert = SAM2PAF(infile, tempfile.name)
        try:
            convert()
            assert False
        except ValueError:
            assert True


def test_conv_extra():
    # calls with SAM/summary/None ar extra_fields argument
    # call with pri_only = False/True
    # Input contains a cigar with all fields MINSH=X 
    infile = bioconvert_data("test_sam2paf_extra.sam", where)

    with TempFile(suffix=".paf") as tempfile:
        convert = SAM2PAF(infile, tempfile.name)
        convert(extra_fields="SAM", pri_only=False)
        convert(extra_fields="SAM", pri_only=True)

        convert(extra_fields="summary")
        convert(extra_fields=None)
