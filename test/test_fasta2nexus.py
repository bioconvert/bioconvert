import pytest
from easydev import TempFile, md5

from bioconvert.fasta2nexus import FASTA2NEXUS

from . import test_dir

@pytest.mark.parametrize("method", FASTA2NEXUS.available_methods)
def test_fa2nx_biopython(method):
    infile = f"{test_dir}/data/fasta/{method}.fasta"
    outfile = f"{test_dir}/data/nexus/{method}.nexus"
    with TempFile(suffix=".nx") as tempfile:
        converter = FASTA2NEXUS(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)