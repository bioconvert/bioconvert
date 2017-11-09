import pytest
from bioconvert.gz2dsrc import GZ2DSRC
from easydev import TempFile, md5


def test_gz2dsrc():
    """
    Test that fasta gz file is converted as expected to a fasta .dsrc file
    """
    from bioconvert import bioconvert_data
    in_gz = bioconvert_data("SP1.fq.gz")
    exp_dsrc = bioconvert_data("SP1.dsrc")
    with TempFile(suffix=".dsrc") as tempfile:
        converter = GZ2DSRC(in_gz, tempfile.name)
        converter()

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(exp_dsrc)
