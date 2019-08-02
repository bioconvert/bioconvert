import pytest
from easydev import TempFile

from bioconvert import bioconvert_data
from bioconvert.cram2fasta import CRAM2FASTA


@pytest.mark.skipif(CRAM2FASTA._method_samtools.is_disabled, reason="missing dependencies")
def test_conv():
    infile = bioconvert_data("test_measles.cram")

    with TempFile(suffix=".fasta") as tempfile:
        convert = CRAM2FASTA(infile, tempfile.name )
        convert(method="samtools")

    for ext in ['gz', 'bz2', 'dsrc']:
        with TempFile(suffix=".fasta.{}".format(ext)) as tempfile:
            convert = CRAM2FASTA(infile, tempfile.name )
            convert(method="samtools")

    infile = bioconvert_data("test_measles_unpaired.sorted.cram")
    with TempFile(suffix=".fasta") as tempfile:
        convert = CRAM2FASTA(infile, tempfile.name )
        convert(method="samtools")

    for ext in ['gz', 'bz2', 'dsrc']:
        with TempFile(suffix=".fasta.{}".format(ext)) as tempfile:
            convert = CRAM2FASTA(infile, tempfile.name )
            convert(method="samtools")
