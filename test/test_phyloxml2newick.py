import pytest
from easydev import TempFile, md5

from bioconvert.phyloxml2newick import PHYLOXML2NEWICK

from . import test_dir


@pytest.mark.parametrize("method", PHYLOXML2NEWICK.available_methods)
def test_xml2nw_biopython(method):
    infile = f"{test_dir}/data/phyloxml/{method}.xml"
    outfile = f"{test_dir}/data/newick/{method}.newick"
    with TempFile(suffix=".newick") as tempfile:
        converter = PHYLOXML2NEWICK(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
