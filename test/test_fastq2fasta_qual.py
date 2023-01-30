from bioconvert.fastq2fasta_qual import FASTQ2FASTA_QUAL
from bioconvert import TempFile, md5

from . import test_dir

# TODO: Add test of the unwrap_fasta method
def test_conv():
    infile = f"{test_dir}/data/gz/measles_R2.fastq.gz"
    with TempFile(suffix=".fasta") as fout1, TempFile(suffix=".qual") as fout2:
        c = FASTQ2FASTA_QUAL(infile, (fout1.name, fout2.name))
        c()


def test_conv():
    infile = f"{test_dir}/data/fastq/ERR.fastq"
    with TempFile(suffix=".fasta") as fout1, TempFile(suffix=".qual") as fout2:
        c = FASTQ2FASTA_QUAL(infile, (fout1.name, fout2.name))
        c()
