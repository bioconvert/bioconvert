import os
from bioconvert.sam2paf import SAM2PAF
from bioconvert import TempFile, md5

from . import test_dir


def test_conv():
    infile = f"{test_dir}/data/sam/test_sam2paf_v1.sam"
    outfile = f"{test_dir}/data/paf/test_sam2paf_v1.paf"
    checksum = md5(outfile)

    with TempFile(suffix=".paf") as tempfile:
        convert = SAM2PAF(infile, tempfile.name)
        convert()
        assert checksum == md5(tempfile.name)
        assert convert.skipped == 17


def test_bad1_input():
    infile = f"{test_dir}/data/sam/test_sam2paf_bad1.sam"
    with TempFile(suffix=".paf") as tempfile:
        convert = SAM2PAF(infile, tempfile.name)
        try:
            convert()
            assert False
        except ValueError:
            assert True


def test_bad2_input():
    infile = f"{test_dir}/data/sam/test_sam2paf_bad2.sam"
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
    infile = f"{test_dir}/data/sam/test_sam2paf_extra.sam"

    with TempFile(suffix=".paf") as tempfile:
        convert = SAM2PAF(infile, tempfile.name)
        convert(extra_fields="SAM", pri_only=False)
        convert(extra_fields="SAM", pri_only=True)

        convert(extra_fields="summary")
        convert(extra_fields=None)
