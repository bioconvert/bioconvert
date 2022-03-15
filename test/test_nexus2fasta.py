import pytest
from easydev import TempFile, md5

from bioconvert.nexus2fasta import NEXUS2FASTA

from . import test_dir

@pytest.mark.parametrize("method", NEXUS2FASTA.available_methods)
def test_nx2fa_biopython(method):
    if method == "goalign":
        infile = f"{test_dir}/data/nexus/goalign.nexus"
        outfile = f"{test_dir}/data/fasta/goalign.fasta"
        with TempFile(suffix=".fasta") as tempfile:
            converter = NEXUS2FASTA(infile, tempfile.name)
            converter(method=method)
            # Check that the output is correct with a checksum
            assert md5(tempfile.name) == md5(outfile)
    else:
        infile = f"{test_dir}/data/nexus/nexus2fasta_{method}.nexus"
        outfile = f"{test_dir}/data/fasta/nexus2fasta_{method}.fasta"
        with TempFile(suffix=".fasta") as tempfile:
            converter = NEXUS2FASTA(infile, tempfile.name)
            converter(method=method)

            # Check that the output is correct with a checksum
            assert md5(tempfile.name) == md5(outfile)