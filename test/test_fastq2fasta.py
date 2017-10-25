from bioconvert.fastq2fasta import Fastq2Fasta
from bioconvert import bioconvert_data
from easydev import TempFile, md5


def test_conv():
    infile = bioconvert_data("test_fastq2fasta_v1.fastq")

    with TempFile(suffix=".fasta") as tempfile:
        convert = Fastq2Fasta(infile, tempfile.name)
        convert()
        print(md5(bioconvert_data("test_fastq2fasta_v1.fasta")))

