from bioconvert.vcf2bed import VCF2BED
from bioconvert import bioconvert_data
from easydev import TempFile, md5

def test_conv():
    infile = bioconvert_data("test_vcf2bcf_v1.vcf")
    outfile = bioconvert_data("test_vcf2bed_v1.bed")
    with TempFile(suffix=".bed") as tempfile:
        convert = VCF2BED(infile, tempfile.name)
        convert(method="awk")
        assert md5(tempfile.name) == md5(outfile)
