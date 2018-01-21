import os
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.nexus2phyloxml import NEXUS2PHYLOXML
import pytest


@pytest.mark.parametrize("method", NEXUS2PHYLOXML.available_methods)
def test_nx2xml_biopython(method):
    infile = bioconvert_data(method+".nexus")
    outfile = bioconvert_data(method+".xml")
    with TempFile(suffix=".phyloxml") as tempfile:
        converter = NEXUS2PHYLOXML(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


