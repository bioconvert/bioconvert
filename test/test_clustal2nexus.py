import os
import pytest
from bioconvert import TempFile, md5

from bioconvert.clustal2nexus import CLUSTAL2NEXUS

from . import test_dir


@pytest.mark.skipif(
    CLUSTAL2NEXUS._method_goalign.is_disabled, reason="missing dependencies"
)
def test_clustal2nexus_goalign():
    infile = f"{test_dir}/data/clustal/goalign.clustal"
    outfile = f"{test_dir}/data/nexus/goalign.nexus"
    with TempFile(suffix=".nexus") as tempfile:
        converter = CLUSTAL2NEXUS(infile, tempfile.name)
        converter(method="goalign")

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
