from bioconvert.fastq2fasta import Fastq2Fasta
from bioconvert import bioconvert_data
from easydev import TempFile, md5
import pytest
from bioconvert.genbank2fasta import GENBANK2FASTA

where = "testing/genbank2fasta"


@pytest.mark.parametrize("method", GENBANK2FASTA.available_methods)
def test_conv(method):
    infile = bioconvert_data("JB409847.gbk", where)

    with TempFile(suffix=".fasta") as tempfile:
        converter = GENBANK2FASTA(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        if method == "biopython":
            assert md5(tempfile.name) == "3c679a17f69cc26afa899b6cf571dd13"
        elif method == "squizz":
            assert md5(tempfile.name) == "40ca920aae3a3481eb2fbe365eadbf4a"
        else:
            raise NotImplementedError


