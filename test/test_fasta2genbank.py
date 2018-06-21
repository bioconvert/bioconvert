from bioconvert import bioconvert_data
from easydev import TempFile, md5
from bioconvert.readers.genbank import Genbank
from bioconvert.readers.fasta import Fasta

import pytest
from bioconvert.fasta2genbank import FASTA2GENBANK



@pytest.mark.parametrize("method", FASTA2GENBANK.available_methods)
def test_conv(method):
    infile = bioconvert_data("test_measles.fa")

    with TempFile(suffix=".gbk") as tempfile:
        converter = FASTA2GENBANK(infile, tempfile.name)
        converter(method=method)

        reader_fasta = Fasta(infile)
        reader_gbk = Genbank(tempfile.name)

        for i in range(10):
        	print(i)

        for fasta_entry, gbk_entry in zip(reader_fasta.read(), reader_gbk.read()):
        	print(fasta_entry)
        	print(gbk_entry)
        	assert fasta_entry["id"] == gbk_entry["LOCUS"]["id"]
        	assert fasta_entry["comment"] in gbk_entry["DEFINITION"]
        	assert fasta_entry["value"] == gbk_entry["ORIGIN"].upper()

        # assert md5(tempfile.name) in out_md5s, "Incorect gbk output for method {}".format(method)

