import os
import pytest
from bioconvert import TempFile, md5

from bioconvert.clustal2fasta import CLUSTAL2FASTA

from . import test_dir


def test_clustal2fasta_biopython():
    infile = f"{test_dir}/data/clustal/biopython.clustal"
    outfile = f"{test_dir}/data/fasta/biopython.fasta"
    with TempFile(suffix=".fasta") as tempfile:
        converter = CLUSTAL2FASTA(infile, tempfile.name)
        converter(method="biopython")

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


@pytest.mark.skipif(
    CLUSTAL2FASTA._method_squizz.is_disabled, reason="missing dependencies"
)
def test_clustal2fasta_squizz():
    infile = f"{test_dir}/data/clustal/squizz.clustal"
    outfile = f"{test_dir}/data/fasta/squizz.fasta"
    with TempFile(suffix=".fasta") as tempfile:
        converter = CLUSTAL2FASTA(infile, tempfile.name)
        converter(method="squizz")

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


@pytest.mark.skipif(
    CLUSTAL2FASTA._method_goalign.is_disabled, reason="missing dependencies"
)
def test_clustal2fasta_goalign():
    infile = f"{test_dir}/data/clustal/goalign.clustal"
    outfile = f"{test_dir}/data/fasta/goalign.fasta"
    with TempFile(suffix=".fasta") as tempfile:
        converter = CLUSTAL2FASTA(infile, tempfile.name)
        converter(method="goalign")

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
