import pytest
from easydev import TempFile, md5

from bioconvert.wig2bed import WIG2BED

from . import test_dir


@pytest.mark.parametrize("method", WIG2BED.available_methods)
def test_nx2xml_biopython(method):
    infile = f"{test_dir}/data/wig/test_wig2bed.wig"
    outfile = f"{test_dir}/data/bed/test_wig2bed.bed"
    with TempFile(suffix=".phyloxml") as tempfile:
        converter = WIG2BED(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
