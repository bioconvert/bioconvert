import pytest
from bioconvert import TempFile

from bioconvert.cram2fastq import CRAM2FASTQ

from . import test_dir


@pytest.mark.skipif(
    CRAM2FASTQ._method_samtools.is_disabled, reason="missing dependencies"
)
def test_conv():
    infile = f"{test_dir}/data/cram/test_measles.cram"

    with TempFile(suffix=".fastq") as tempfile:
        convert = CRAM2FASTQ(infile, tempfile.name)
        convert(method="samtools")

    for ext in ["gz", "bz2", "dsrc"]:
        with TempFile(suffix=".fastq.{}".format(ext)) as tempfile:
            convert = CRAM2FASTQ(infile, tempfile.name)
            convert(method="samtools")

    infile = f"{test_dir}/data/cram/test_measles_unpaired.sorted.cram"
    with TempFile(suffix=".fastq") as tempfile:
        convert = CRAM2FASTQ(infile, tempfile.name)
        convert(method="samtools")

    for ext in ["gz", "bz2", "dsrc"]:
        with TempFile(suffix=".fastq.{}".format(ext)) as tempfile:
            convert = CRAM2FASTQ(infile, tempfile.name)
            convert(method="samtools")
