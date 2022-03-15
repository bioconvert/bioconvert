import os
import pytest
from easydev import TempFile, md5

from bioconvert.clustal2phylip import CLUSTAL2PHYLIP

from . import test_dir

skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ
                                and os.environ['TRAVIS_PYTHON_VERSION'].startswith("2"), reason="On travis")


@skiptravis
def test_clustal2phylip_biopython():
    infile = f"{test_dir}/data/clustal/biopython.clustal"
    outfile = f"{test_dir}/data/phylip/biopython.phylip"
    with TempFile(suffix=".phylip") as tempfile:
        converter = CLUSTAL2PHYLIP(infile, tempfile.name)
        converter(method='biopython')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)

@skiptravis
@pytest.mark.skipif(CLUSTAL2PHYLIP._method_squizz.is_disabled, reason="missing dependencies")
def test_clustal2phylip_squizz():
    infile = f"{test_dir}/data/clustal/squizz.clustal"
    outfile = f"{test_dir}/data/phylip/squizz.phylip"
    with TempFile(suffix=".phylip") as tempfile:
        converter = CLUSTAL2PHYLIP(infile, tempfile.name)
        converter(method='squizz')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)