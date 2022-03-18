import os
import pytest
import bioconvert
from bioconvert.core.utils import generate_outfile_name

from . import test_dir


def test_generate_outfile_name():
    assert generate_outfile_name("foo.fasta", "phylip") == "foo.phylip"
    assert (
        generate_outfile_name("/foo/bar.ext.fasta", "clustal") == "/foo/bar.ext.clustal"
    )
