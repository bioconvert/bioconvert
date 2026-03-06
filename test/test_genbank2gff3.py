from bioconvert import TempFile, md5
import pytest
from bioconvert.genbank2gff3 import GENBANK2GFF3
from bioconvert.core.registry import Registry

from . import test_dir

# This test fails with wierd pytest error message related
# to io module ?
# FIXME: https://github.com/bioconvert/bioconvert/issues/143
def test_conv():
    infile = f"{test_dir}/data/genbank/biocode.gb"
    outfile = f"{test_dir}/data/GFF3/biocode.gff"

    with TempFile(suffix=".gff") as tempfile:
        converter = GENBANK2GFF3(infile, tempfile.name)
        converter(method="biocode")
        assert md5(tempfile.name) == md5(outfile)


def test_gbff_extension_recognized():
    """Test that .gbff extension is recognized as a valid input format for genbank2gff3."""
    registry = Registry()
    assert (("gbff",), ("gff3",)) in registry._ext_registry, (
        ".gbff extension should be recognized for genbank2gff3 conversion"
    )

