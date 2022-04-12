import os
import pytest
from easydev import TempFile, md5

from bioconvert.bigbed2bed import BIGBED2BED
import pytest

from . import test_dir


@pytest.mark.parametrize("method", BIGBED2BED.available_methods)
def test_bigwig2bedgraph_ucsc(method):
    infile = f"{test_dir}/data/bigbed/test_pybigwig.bigbed"
    outfile = f"{test_dir}/data/bed/test_pybigwig.bed"
    with TempFile(suffix=".bed") as tempfile:
        converter = BIGBED2BED(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
