from bioconvert.io.maf import MAF
from bioconvert.io import maf
from easydev import TempFile

from bioconvert import bioconvert_data


filename = bioconvert_data("test_maf2sam.maf")


def test_read_maf():

    with TempFile(suffix=".sam") as fout:
        maf = MAF(filename, fout.name)
        #maf.count_insertions()
        maf.to_sam()
        
def test_others():
    maf.mapqFromProb("0.5")
