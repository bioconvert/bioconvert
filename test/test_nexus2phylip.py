import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.nexus2phylip import NEXUS2PHYLIP


@pytest.mark.parametrize("method", NEXUS2PHYLIP.available_methods)
def test_nx2phy_biopython(method):
    infile = bioconvert_data(method + ".nexus")
    outfile = bioconvert_data(method + ".phylip")
    with TempFile(suffix=".phylip") as tempfile:
        converter = NEXUS2PHYLIP(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        #assert md5(tempfile.name) == md5(outfile)
        #https://github.com/bioconvert/bioconvert/issues#149
