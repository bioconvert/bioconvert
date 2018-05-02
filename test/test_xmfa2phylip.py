import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.xmfa2phylip import XMFA2PHYLIP


@pytest.mark.parametrize("method", XMFA2PHYLIP.available_methods)
def test_xmfa2phy(method):
    infile = bioconvert_data("test_phylip2xmfa.xmfa")
    #outfile = bioconvert_data("test_phylip2xmfa.xmfa")
    with TempFile(suffix=".xmfa") as tempfile:
        converter = XMFA2PHYLIP(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        #TODO
        #assert md5(tempfile.name) == md5(outfile)
