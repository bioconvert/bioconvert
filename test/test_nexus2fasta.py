import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.nexus2fasta import NEXUS2FASTA


@pytest.mark.parametrize("method", NEXUS2FASTA.available_methods)
def test_nx2fa_biopython(method):
    if method == "goalign":
        infile = bioconvert_data("goalign.nexus")
        outfile = bioconvert_data("goalign.fasta")
        with TempFile(suffix=".fasta") as tempfile:
            converter = NEXUS2FASTA(infile, tempfile.name)
            converter(method=method)
            # Check that the output is correct with a checksum
            assert md5(tempfile.name) == md5(outfile)
    else:
        infile = bioconvert_data("nexus2fasta_"+ method + ".nexus")
        outfile = bioconvert_data("nexus2fasta_" + method + ".fasta")
        with TempFile(suffix=".fasta") as tempfile:
            converter = NEXUS2FASTA(infile, tempfile.name)
            converter(method=method)

            # Check that the output is correct with a checksum
            assert md5(tempfile.name) == md5(outfile)
