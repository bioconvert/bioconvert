import pytest
from easydev import TempFile

from bioconvert.cram2fasta import CRAM2FASTA

from . import test_dir


@pytest.mark.skipif(
    CRAM2FASTA._method_samtools.is_disabled, reason="missing dependencies"
)
def test_conv():
    infile = f"{test_dir}/data/cram/test_measles.cram"

    with TempFile(suffix=".fasta") as tempfile:
        convert = CRAM2FASTA(infile, tempfile.name)
        convert(method="samtools")

    for ext in ["gz", "bz2", "dsrc"]:
        with TempFile(suffix=".fasta.{}".format(ext)) as tempfile:
            convert = CRAM2FASTA(infile, tempfile.name)
            convert(method="samtools")

    infile = f"{test_dir}/data/cram/test_measles_unpaired.sorted.cram"
    with TempFile(suffix=".fasta") as tempfile:
        convert = CRAM2FASTA(infile, tempfile.name)
        convert(method="samtools")

    for ext in ["gz", "bz2", "dsrc"]:
        with TempFile(suffix=".fasta.{}".format(ext)) as tempfile:
            convert = CRAM2FASTA(infile, tempfile.name)
            convert(method="samtools")
