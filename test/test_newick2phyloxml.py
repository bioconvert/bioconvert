import os
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.newick2phyloxml import NEWICK2PHYLOXML
import pytest


@pytest.mark.parametrize("method", NEWICK2PHYLOXML.available_methods)
def test_nw2xml_biopython(method):
    infile = bioconvert_data(method+".newick")
    outfile = bioconvert_data(method+".xml")
    with TempFile(suffix=".xml") as tempfile:
        converter = NEWICK2PHYLOXML(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
