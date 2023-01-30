import pytest
from bioconvert import TempFile, md5

from bioconvert.newick2phyloxml import NEWICK2PHYLOXML

from . import test_dir


@pytest.mark.parametrize("method", NEWICK2PHYLOXML.available_methods)
def test_nw2xml_biopython(method):
    infile = f"{test_dir}/data/newick/{method}.newick"
    outfile = f"{test_dir}/data/phyloxml/{method}.xml"
    with TempFile(suffix=".xml") as tempfile:
        converter = NEWICK2PHYLOXML(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
