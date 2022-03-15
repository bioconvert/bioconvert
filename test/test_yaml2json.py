from bioconvert.yaml2json import YAML2JSON
from easydev import TempFile, md5
import pytest
import os

from . import test_dir

skiptravis = pytest.mark.skipif( "TRAVIS_PYTHON_VERSION" in os.environ and
    os.environ['TRAVIS_PYTHON_VERSION'].startswith("2"), reason="On travis")


@skiptravis
def test_conv():
    infile = f"{test_dir}/data/yaml/test_v1.yaml"
    expected_outile = f"{test_dir}/data/json/test_v1.json"
    with TempFile(suffix=".json") as tempfile:
        convert = YAML2JSON(infile, tempfile.name)
        convert()

        # Check that the output is correct with a checksum
        # Note that we cannot test the md5 on a gzip file but only 
        # on the original data. This check sum was computed
        # fro the unzipped version of bioconvert/data/measles.bed
        assert md5(tempfile.name) == md5(expected_outile)