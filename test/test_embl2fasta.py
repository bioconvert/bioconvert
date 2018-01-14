from bioconvert import bioconvert_data
from easydev import TempFile, md5
import pytest
from bioconvert.embl2fasta import EMBL2FASTA



@pytest.mark.parametrize("method", EMBL2FASTA.available_methods)
def test_conv(method):
    infile = bioconvert_data("JB409847.embl")

    with TempFile(suffix=".gbk") as tempfile:
        converter = EMBL2FASTA(infile, tempfile.name)
        converter(method=method)

        # FIXME
        # need to check md5


