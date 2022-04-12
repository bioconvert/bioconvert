import os
import pytest
from easydev import TempFile, md5

from bioconvert.bedgraph2bigwig import BEDGRAPH2BIGWIG
import pytest

from . import test_dir

url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"


@pytest.mark.parametrize("method", BEDGRAPH2BIGWIG.available_methods)
def test_bigwig2bedgraph_ucsc(method):
    infile = f"{test_dir}/data/bedgraph/ucsc.bedgraph"
    outfile = f"{test_dir}/data/bigwig/ucsc.bigwig"
    with TempFile(suffix=".bigwig") as tempfile:
        converter = BEDGRAPH2BIGWIG(infile, tempfile.name)
        converter(method=method, chrom_sizes=url)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
