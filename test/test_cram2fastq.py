import pytest
from easydev import TempFile

from bioconvert import bioconvert_data
from bioconvert.cram2fastq import CRAM2FASTQ


@pytest.mark.skipif(CRAM2FASTQ._method_samtools.is_disabled, reason="missing dependencies")
def test_conv():
    infile = bioconvert_data("test_measles.cram")
    outfile = bioconvert_data("test_measles.sam")

    with TempFile(suffix=".fastq") as tempfile:
        convert = CRAM2FASTQ(infile, tempfile.name )
        convert(method="samtools")
