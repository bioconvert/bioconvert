import os
import subprocess

import pytest
from easydev import TempFile
from easydev import md5

from bioconvert.gz2dsrc import GZ2DSRC

from . import test_dir


#@pytest.mark.skipif(len(GZ2DSRC.available_methods) == 0, reason="missing dependencies")
def test_gz2dsrc():
    """
    Test that fastq gz file is converted as expected to a fastq .dsrc file
    """
    in_gz = f"{test_dir}/data/gz/test_SP1.fq.gz"
    exp_fq = f"{test_dir}/data/fastq/exp_SP1.fq"
    with TempFile(suffix=".dsrc") as tempfile:
        converter = GZ2DSRC(in_gz, tempfile.name)
        converter()

        # uncompress the created dsrc file, and compare uncompressed file
        # to the expected one. We do not directly compare dsrc files as
        # it depends on the dsrc version used...
        assert os.path.isfile(tempfile.name)
        tmp_fq = tempfile.name + ".fq"
        cmd = "dsrc d {} {}".format(tempfile.name, tmp_fq)
        subprocess.call(cmd.split())

        # Check that the output is correct with a checksum
        assert md5(tmp_fq) == md5(exp_fq)
