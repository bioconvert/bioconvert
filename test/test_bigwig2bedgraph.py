import os
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.bigwig2bedgraph import BIGWIG2BEDGRAPH
import pytest


skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ
                                and os.environ['TRAVIS_PYTHON_VERSION'].startswith("2"), reason="On travis")


@pytest.mark.parametrize("method", BIGWIG2BEDGRAPH.available_methods)
def test_bigwig2bedgraph_ucsc(method):
    infile = bioconvert_data("ucsc.bigwig")
    outfile = bioconvert_data("ucsc.bedgraph")
    with TempFile(suffix=".bedgraph") as tempfile:
        converter = BIGWIG2BEDGRAPH(infile, tempfile.name)
        converter(method='ucsc')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
