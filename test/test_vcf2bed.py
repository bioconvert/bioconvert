from bioconvert.vcf2bed import VCF2BED
from bioconvert import TempFile, md5

from . import test_dir


def test_conv():
    infile = f"{test_dir}/data/vcf/test_vcf2bcf_v1.vcf"
    outfile = f"{test_dir}/data/bed/test_vcf2bed_v1.bed"
    with TempFile(suffix=".bed") as tempfile:
        convert = VCF2BED(infile, tempfile.name)
        convert(method="awk")
        assert md5(tempfile.name) == md5(outfile)
