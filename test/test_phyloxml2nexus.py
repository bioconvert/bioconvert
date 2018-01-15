import os
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.phyloxml2nexus import PHYLOXML2NEXUS
import pytest


@pytest.mark.parametrize("method", PHYLOXML2NEXUS.available_methods)
def test_xml2nx_biopython(method):
    infile = bioconvert_data(method+".xml")
    outfile = bioconvert_data(method+".nexus")
    with TempFile(suffix=".nexus") as tempfile:
        converter = PHYLOXML2NEXUS(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
