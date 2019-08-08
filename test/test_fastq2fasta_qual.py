from bioconvert.fastq2fasta_qual import FASTQ2FASTA_QUAL
from bioconvert import bioconvert_data
from easydev import TempFile, md5


# TODO: Add test of the unwrap_fasta method
def test_conv():
    infile = bioconvert_data("test_fastq2fasta_v1.fastq")

    expected_outfile = bioconvert_data("test_fastq2fasta_v1.fasta")
    with TempFile(suffix=".fasta") as fout1, TempFile(suffix=".qual") as fout2:
        c = FASTQ2FASTA_QUAL(infile, (fout1.name, fout2.name))
        c()


