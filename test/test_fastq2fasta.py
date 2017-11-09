from bioconvert.fastq2fasta import Fastq2Fasta
from bioconvert import bioconvert_data
from easydev import TempFile, md5


# TODO: Add test of the unwrap_fasta method
def test_conv():
    infile = bioconvert_data("test_fastq2fasta_v1.fastq")

    # this md5sum contains @ instead of > for the identifiers
    # outfile = bioconvert_data("test_fastq2fasta_v1.fasta")
    # so for now, we use the md5sum (that would use > signs)
    md5out = "27e2e0ae66fe2e2fb17c6715cafffdf0"

    # One temporary file for the fasta created using the method
    # and one for an unwrapped version.
    # Some methods may output multi-line fasta, so we need to
    # compare md5 sums of unwrapped versions.
    for method in Fastq2Fasta.available_methods:
        with TempFile(suffix=".fasta") as outfile, \
                TempFile(suffix=".fasta") as unwrapped:
            convert = Fastq2Fasta(infile, outfile.name)
            convert()
            Fastq2Fasta.unwrap_fasta(outfile.name, unwrapped.name)
            assert md5(unwrapped.name) == md5out, \
                "{} failed".format(method)


