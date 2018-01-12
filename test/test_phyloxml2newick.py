import os
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.phyloxml2newick import PHYLOXML2NEWICK
import pytest


skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ
                                and os.environ['TRAVIS_PYTHON_VERSION'].startswith("2"), reason="On travis")


@skiptravis
@pytest.mark.parametrize("method", PHYLOXML2NEWICK.available_methods)
def test_xml2nw_biopython(method):
    infile = bioconvert_data("gotree.xml")
    outfile = bioconvert_data("gotree.nw")
    with TempFile(suffix=".nw") as tempfile:
        converter = PHYLOXML2NEWICK(infile, tempfile.name)
        converter(method='gotree')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
