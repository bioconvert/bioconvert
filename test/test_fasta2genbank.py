from easydev import TempFile, md5
from bioconvert.io.genbank import Genbank
from bioconvert.io.fasta import Fasta

import pytest
from bioconvert.fasta2genbank import FASTA2GENBANK

from . import test_dir


@pytest.mark.parametrize("method", FASTA2GENBANK.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/fasta/test_measles.fa"

    with TempFile(suffix=".gbk") as tempfile:
        converter = FASTA2GENBANK(infile, tempfile.name)
        converter(method=method)

        reader_fasta = Fasta(infile)
        reader_gbk = Genbank(tempfile.name)

        for fasta_entry, gbk_entry in zip(reader_fasta.read(), reader_gbk.read()):
            assert fasta_entry["id"] == gbk_entry["LOCUS"]["id"]
            assert fasta_entry["comment"] in gbk_entry["DEFINITION"]
            assert fasta_entry["value"] == gbk_entry["ORIGIN"].upper()
