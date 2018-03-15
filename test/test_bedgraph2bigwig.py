import os
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.bedgraph2bigwig import BEDGRAPH2BIGWIG
import pytest


@pytest.mark.parametrize("method", BEDGRAPH2BIGWIG.available_methods)
def test_bigwig2bedgraph_ucsc(method):
    infile = bioconvert_data("ucsc.bedgraph")
    outfile = bioconvert_data("ucsc.bigwig")
    with TempFile(suffix=".bigwig") as tempfile:
        converter = BEDGRAPH2BIGWIG(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
