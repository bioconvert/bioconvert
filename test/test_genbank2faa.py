import pytest
from bioconvert import TempFile, md5
from bioconvert.genbank2faa import GENBANK2FAA

from . import test_dir


@pytest.mark.parametrize("method", GENBANK2FAA.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/genbank/test_genbank2faa.gbk"
    expected_outfile = f"{test_dir}/data/faa/test_genbank2faa.faa"

    with TempFile(suffix=".faa") as tempfile:
        converter = GENBANK2FAA(infile, tempfile.name)
        converter(method=method)
        assert md5(tempfile.name) == md5(expected_outfile)
