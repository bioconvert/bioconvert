import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.phylip2xmfa import PHYLIP2XMFA


@pytest.mark.parametrize("method", PHYLIP2XMFA.available_methods)
def test_phy2xmfa(method):
    infile = bioconvert_data(method + ".phylip")
    outfile = bioconvert_data("test_phylip2xmfa.xmfa")
    with TempFile(suffix=".xmfa") as tempfile:
        converter = PHYLIP2XMFA(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        #TODO
        #assert md5(tempfile.name) == md5(outfile)
