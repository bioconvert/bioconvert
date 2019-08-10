
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.gff22gff3 import GFF22GFF3


@pytest.mark.parametrize("method", GFF22GFF3.available_methods)
def test_conv(method):
    infile = bioconvert_data("GFF2/gff2_example.gff")
    with TempFile(suffix=".tsv") as tempfile:
        convert = GFF22GFF3(infile, tempfile.name)
        convert(method=method)
