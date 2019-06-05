from bioconvert.fastq2fasta import FASTQ2FASTA
from bioconvert import bioconvert_data
from easydev import TempFile, md5
import pytest
from bioconvert.genbank2fasta import GENBANK2FASTA

from bioconvert.readers.fasta import Fasta
from bioconvert.readers.genbank import Genbank

where = "testing/genbank2fasta"


@pytest.mark.parametrize("method", GENBANK2FASTA.available_methods)
def test_conv(method):
    infile = bioconvert_data("JB409847.gbk", where)

    with TempFile(suffix=".fasta") as tempfile:
        converter = GENBANK2FASTA(infile, tempfile.name)
        converter(method=method)

        # Load both files
        reader_fasta = Fasta(tempfile.name)
        reader_gen = Genbank(infile)

        for entry_fa, entry_gb in zip(reader_fasta.read(), reader_gen.read()):
            assert entry_fa["id"].startswith(entry_gb["LOCUS"]["id"])
            assert entry_fa["comment"] in entry_gb["DEFINITION"]
            assert entry_fa["value"].lower() == entry_gb["ORIGIN"]


