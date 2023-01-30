import os
import pytest
import hashlib
from bioconvert import TempFile, md5

from bioconvert.fasta2clustal import FASTA2CLUSTAL

from . import test_dir


def test_fasta2clustal_biopython():
    infile = f"{test_dir}/data/fasta/biopython.fasta"
    outfile = f"{test_dir}/data/clustal/biopython.clustal"
    with TempFile(suffix=".clustal") as tempfile:
        converter = FASTA2CLUSTAL(infile, tempfile.name)
        converter(method="biopython")

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


@pytest.mark.skipif(
    FASTA2CLUSTAL._method_squizz.is_disabled, reason="missing dependencies"
)
def test_fasta2clustal_squizz():
    infile = f"{test_dir}/data/fasta/squizz.fasta"
    outfile = f"{test_dir}/data/clustal/squizz.clustal"
    with TempFile(suffix=".clustal") as tempfile:
        converter = FASTA2CLUSTAL(infile, tempfile.name)
        converter(method="squizz")

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


@pytest.mark.skipif(
    FASTA2CLUSTAL._method_goalign.is_disabled, reason="missing dependencies"
)
def test_fasta2clustal_goalign():
    infile = f"{test_dir}/data/fasta/goalign.fasta"
    outfile = f"{test_dir}/data/clustal/goalign.clustal"
    with TempFile(suffix=".clustal") as tempfile:
        converter = FASTA2CLUSTAL(infile, tempfile.name)
        converter(method="goalign")

        ## We remove goalign version from the first line
        out = ""
        with open(tempfile.name) as f:
            lines = f.readlines()
            if len(lines) > 0:
                clustal = lines[0].split(" ")
                if len(clustal) > 0:
                    lines[0] = clustal[0] + "\n"
            out = "".join(lines)

        # Check that the output is correct with a checksum
        assert hashlib.md5(out.encode("utf-8")).hexdigest() == md5(outfile)
