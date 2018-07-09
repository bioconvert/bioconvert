from bioconvert import bioconvert_data
from easydev import TempFile, md5
import pytest
from bioconvert.genbank2gff3 import GENBANK2GFF3


"""
@pytest.mark.parametrize("method", GENBANK2GFF.available_methods)
def test_conv(method):
    infile = bioconvert_data(method + ".gb")
    outfile = bioconvert_data(method + ".gff")

    with TempFile(suffix=".gff") as tempfile:
        converter = GENBANK2GFF(infile, tempfile.name)
        converter(method=method)
        assert md5(tempfile.name) == md5(outfile)
"""
