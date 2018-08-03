import os
import subprocess
import hashlib
import pytest
from easydev import TempFile

from bioconvert.dsrc2gz import DSRC2GZ


@pytest.mark.skipif(len(DSRC2GZ.available_methods) == 0, reason="missing dependencies")
def test_gz2dsrc():
    """
    Test that fastq gz file is converted as expected to a fastq .dsrc file
    """
    from bioconvert import bioconvert_data
    infile = bioconvert_data("test_SP1.fq.dsrc")

    with TempFile(suffix=".fq.gz") as tempfile:
        converter = DSRC2GZ(infile, tempfile.name)
        converter()

        # uncompress the createdfile, and compare uncompressed file
        # to the expected md5. We do not directly compare dsrc or gz files as
        # it is not deterministic
        assert os.path.isfile(tempfile.name)

        cmd = "gunzip -c {}".format(tempfile.name)
        res = subprocess.check_output(cmd, shell=True)

        # Check that the output is correct with a checksum
        assert hashlib.md5(res).hexdigest() == "d41d8cd98f00b204e9800998ecf8427e"
