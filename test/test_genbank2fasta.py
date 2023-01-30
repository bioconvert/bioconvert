from bioconvert.fastq2fasta import FASTQ2FASTA
from bioconvert import TempFile, md5
import pytest
from bioconvert.genbank2fasta import GENBANK2FASTA

from bioconvert.io.fasta import Fasta
from bioconvert.io.genbank import Genbank

from . import test_dir


@pytest.mark.parametrize("method", GENBANK2FASTA.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/genbank/JB409847.gbk"

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
