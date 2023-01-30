import pytest
from bioconvert import TempFile, md5

from bioconvert.gff32gff2 import GFF32GFF2

from . import test_dir


@pytest.mark.parametrize("method", GFF32GFF2.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/GFF3/gff3_example.gff"
    with TempFile(suffix=".tsv") as tempfile:
        convert = GFF32GFF2(infile, tempfile.name)
        convert(method=method)
