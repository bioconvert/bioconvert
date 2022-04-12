from bioconvert.maf2sam import MAF2SAM
from easydev import TempFile, md5

from . import test_dir


def test_conv():
    infile = f"{test_dir}/data/maf/test_maf2sam.maf"
    outfile = f"{test_dir}/data/sam/test_maf2sam.sam"
    with TempFile(suffix=".sam") as tempfile:
        convert = MAF2SAM(infile, tempfile.name)
        convert(method="python")

        # In the SAM, the version may be different when using other bioconvert
        # version, so we need to get rid of the line that contains the version
        # and program
        data1 = open(outfile).readlines()
        data1 = [x for x in data1 if "bioconvert" not in x]
        data1 = "\n".join(data1)

        data2 = open(tempfile.name).readlines()
        data2 = [x for x in data2 if "bioconvert" not in x]
        data2 = "\n".join(data2)

        assert data1 == data2
