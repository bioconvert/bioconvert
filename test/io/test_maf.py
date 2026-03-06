from bioconvert.io.maf import MAF
from bioconvert.io import maf
from bioconvert import TempFile

from .. import test_dir

filename = f"{test_dir}/data/maf/test_maf2sam.maf"
filename_comments = f"{test_dir}/data/maf/test_maf2sam_comments.maf"


def test_read_maf():

    with TempFile(suffix=".sam") as fout:
        maf = MAF(filename, fout.name)
        # maf.count_insertions()
        maf.to_sam()


def test_read_maf_with_comments():
    """Regression test: MAF files with comment lines (#) should not fail."""
    with TempFile(suffix=".sam") as fout:
        m = MAF(filename_comments, fout.name)
        m.to_sam()


def test_others():
    maf.mapqFromProb("0.5")
