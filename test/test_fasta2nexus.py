import os
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.fasta2nexus import FASTA2NEXUS
import pytest


skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ
                                and os.environ['TRAVIS_PYTHON_VERSION'].startswith("2"), reason="On travis")


@skiptravis
@pytest.mark.parametrize("method", FASTA2NEXUS.available_methods)
def test_fa2nx_biopython(method):
    infile = bioconvert_data("goalign.fa")
    outfile = bioconvert_data("goalign.nx")
    with TempFile(suffix=".nx") as tempfile:
        converter = FASTA2NEXUS(infile, tempfile.name)
        converter(method='goalign')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


