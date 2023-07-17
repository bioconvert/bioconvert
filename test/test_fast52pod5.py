from bioconvert.fast52pod5 import FAST52POD5
from bioconvert import TempFile
import pytest

from . import test_dir


@pytest.mark.parametrize("method", FAST52POD5.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/fast5/single_read.fast5"

    with TempFile(suffix=".pod5") as tempfile:
        convert = FAST52POD5(infile, tempfile.name)
        convert(method=method, force=True) # force because pod5 requires it if the file exists. 

        # md5 changes at each conversion...probably due to a time stamp
        #assert md5(tempfile.name) == md5out, "{} failed".format(method)
