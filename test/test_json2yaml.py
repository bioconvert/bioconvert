import pytest
import os

from bioconvert.json2yaml import JSON2YAML
from easydev import TempFile, md5

from . import test_dir


def test_conv():
    infile = f"{test_dir}/data/json/test_v1.json"
    expected_outile = f"{test_dir}/data/yaml/test_v1_nocomments.yaml"
    with TempFile(suffix=".yaml") as tempfile:
        convert = JSON2YAML(infile, tempfile.name)
        convert()

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(expected_outile)
