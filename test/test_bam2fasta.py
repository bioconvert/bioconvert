import pytest

from bioconvert.bam2fasta import BAM2FASTA
from bioconvert import bioconvert_data
from easydev import TempFile, md5


#@pytest.mark.skipif(BAM2FASTA._method_sammtools.is_disabled, reason="missing dependencies")
#def test_method_bamtools():
#    infile = bioconvert_data("test_measles.sorted.bam")
#    with TempFile(suffix=".fa") as tempfile:
#        convert = BAM2FASTA(infile, tempfile.name)
#        convert(method="bamtools")
#
#        # Check that the output is correct with a checksum
#        # Note that we cannot test the md5 on a gzip file but only 
#        # on the original data. This check sum was computed
#        # fro the unzipped version of biokit/data/converters/measles.bed
#        assert md5(tempfile.name) == "ea5511c3c8913626be152609887c8c4d"


@pytest.mark.skipif(BAM2FASTA._method_samtools.is_disabled, reason="missing dependencies")
def test_method_samtools():
    infile = bioconvert_data("test_measles.sorted.bam")

    with TempFile(suffix=".fa") as tempfile:
        convert = BAM2FASTA(infile, tempfile.name)
        convert(method="samtools")
        # samtools 1.6 / hstlib 1.6 gives different results on travis and
        # locally
        assert md5(tempfile.name.replace(".", "_1.")) in [
            "9242d127969a089ddeedbc2002c62686"]
        assert md5(tempfile.name.replace(".", "_2.")) in [
            "b753ad368c9614130884acb29861bd23"]

    # Test compression 
    with TempFile(suffix=".fasta.gz") as tempfile:
        convert = BAM2FASTA(infile, tempfile.name)
        convert(method="samtools")
        # no check, just running 

    # Test compression 
    with TempFile(suffix=".fasta.bz2") as tempfile:
        convert = BAM2FASTA(infile, tempfile.name)
        convert(method="samtools")
        # no check, just running 
