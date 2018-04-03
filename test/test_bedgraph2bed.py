import os
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.bedgraph2bed import BEDGRAPH2BED
import pytest


@pytest.mark.parametrize("method", BEDGRAPH2BED.available_methods)
def test_bedgraph2bed(method):
    infile = bioconvert_data("test_bedgraph2bed.bedgraph")
    with TempFile(suffix=".bigwig") as tempfile:
        converter = BEDGRAPH2BED(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == "a8cc8b0fd2f2fd028424dc8969a0b8b6"
