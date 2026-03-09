import pytest

from bioconvert.fast52fastq import FAST52FASTQ
from bioconvert import TempFile

from . import test_dir


@pytest.mark.parametrize("method", FAST52FASTQ.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/fast5/basecalled.fast5"

    with TempFile(suffix=".fastq") as tempfile:
        convert = FAST52FASTQ(infile, tempfile.name)
        convert(method=method)

        # Verify output is a valid FASTQ file
        with open(tempfile.name, "r") as fq:
            lines = fq.read().splitlines()

        assert len(lines) >= 4, "FASTQ output should have at least 4 lines per read"
        assert lines[0].startswith("@"), "First line should start with '@'"
        assert lines[2] == "+", "Third line should be '+'"
        assert len(lines[1]) == len(lines[3]), "Sequence and quality lengths should match"


@pytest.mark.parametrize("method", FAST52FASTQ.available_methods)
def test_conv_no_basecall(method):
    """Test that a ValueError is raised when the FAST5 has no basecalled data."""
    infile = f"{test_dir}/data/fast5/single_read.fast5"

    with TempFile(suffix=".fastq") as tempfile:
        convert = FAST52FASTQ(infile, tempfile.name)
        with pytest.raises(ValueError, match="No basecalled FASTQ data found"):
            convert(method=method)
