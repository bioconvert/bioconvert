import os
import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.fasta2twobit import FASTA2TWOBIT
import pytest


skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ
                                and os.environ['TRAVIS_PYTHON_VERSION'].startswith("2"), reason="On travis")


@pytest.mark.parametrize("method", FASTA2TWOBIT.available_methods)
def test_fasta2twobit_ucsc(method):
    infile = bioconvert_data("ucsc.fasta")
    outfile = bioconvert_data("ucsc.2bit")
    with TempFile(suffix=".2bit") as tempfile:
        converter = FASTA2TWOBIT(infile, tempfile.name)
        converter(method='ucsc')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
