import os
import pytest
from easydev import TempFile, md5

from bioconvert.bigwig2bedgraph import BIGWIG2BEDGRAPH
import pytest

from . import test_dir

@pytest.mark.parametrize("method", BIGWIG2BEDGRAPH.available_methods)
def test_bigwig2bedgraph_ucsc(method):
    infile = f"{test_dir}/data/bigwig/ucsc.bigwig"
    outfile = f"{test_dir}/data/bedgraph/ucsc.bedgraph"
    with TempFile(suffix=".bedgraph") as tempfile:
        converter = BIGWIG2BEDGRAPH(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)