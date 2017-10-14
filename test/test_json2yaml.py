import pytest
import os

from bioconvert.json2yaml import JSON2YAML
from bioconvert import bioconvert_data
from easydev import TempFile, md5


skiptravis = pytest.mark.skipif( "TRAVIS_PYTHON_VERSION" in os.environ and 
    os.environ['TRAVIS_PYTHON_VERSION'].startswith("2"), reason="On travis")


@skiptravis
def test_conv():
    infile = bioconvert_data("test_v1.json")
    expected_outile = bioconvert_data("test_v1_nocomments.yaml")
    with TempFile(suffix=".yaml") as tempfile:
        convert = JSON2YAML(infile, tempfile.name)
        convert()

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(expected_outile)
