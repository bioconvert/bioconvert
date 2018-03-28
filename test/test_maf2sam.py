from bioconvert.maf2sam import MAF2SAM
from bioconvert import bioconvert_data
from easydev import TempFile, md5


def test_conv():
    infile = bioconvert_data("test_maf2sam.maf")
    outfile = bioconvert_data("test_maf2sam.sam")
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
        data2  = [x for x in data2 if "bioconvert" not in x]
        data2 = "\n".join(data2)

        assert data1 == data2

        
