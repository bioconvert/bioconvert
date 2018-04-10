from bioconvert import bioconvert_data
from easydev import TempFile, md5

import pytest
from bioconvert.fasta2genbank import FASTA2GENBANK



@pytest.mark.parametrize("method", FASTA2GENBANK.available_methods)
def test_conv(method):
    infile = bioconvert_data("test_measles.fa")
    out_cmp = bioconvert_data("test_measles.gbk")

    with TempFile(suffix=".gbk") as tempfile:
        converter = FASTA2GENBANK(infile, tempfile.name)
        converter(method=method)

        assert md5(tempfile.name) == md5(out_cmp), "Incorect gbk output for method {}".format(method)
        # FIXME check the md5sum



