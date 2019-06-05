import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.gff32gff2 import GFF32GFF2


@pytest.mark.parametrize("method", GFF32GFF2.available_methods)
def test_conv(method):
    infile = bioconvert_data("GFF3/gff3_example.gff")
    with TempFile(suffix=".tsv") as tempfile:
        convert = GFF32GFF2(infile, tempfile.name)
        convert(method=method)
