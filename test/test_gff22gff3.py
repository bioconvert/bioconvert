import pytest
from easydev import TempFile, md5

from bioconvert.gff22gff3 import GFF22GFF3

from . import test_dir


@pytest.mark.parametrize("method", GFF22GFF3.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/GFF2/gff2_example.gff"
    with TempFile(suffix=".tsv") as tempfile:
        convert = GFF22GFF3(infile, tempfile.name)
        convert(method=method)
