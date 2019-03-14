from bioconvert.fastq2qual import FASTQ2QUAL
from bioconvert import bioconvert_data
from bioconvert.core.decorators import make_in_gz_tester, requires
from easydev import TempFile, md5
import pytest


@pytest.mark.parametrize("method", FASTQ2QUAL.available_methods)
def test_conv(method):
    infile = bioconvert_data("test_fastq2fasta_v1.fastq")

    expected_outfile = bioconvert_data("test_fastq2qual_v1.qual")
    with TempFile(suffix=".fasta") as fout:
        FASTQ2QUAL(infile, fout.name)
        # TODO: check md5



