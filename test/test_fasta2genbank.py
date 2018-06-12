from bioconvert import bioconvert_data
from easydev import TempFile, md5

import pytest
from bioconvert.fasta2genbank import FASTA2GENBANK



@pytest.mark.parametrize("method", FASTA2GENBANK.available_methods)
def test_conv(method):
    infile = bioconvert_data("test_measles.fa")
    out_cmp = bioconvert_data("test_measles.gbk")
    out_cmp_biopy = bioconvert_data("test_measles_biopython.gbk")

    out_md5s = [md5(out_cmp), md5(out_cmp_biopy)]



    with TempFile(suffix=".gbk") as tempfile:
        converter = FASTA2GENBANK(infile, tempfile.name)
        converter(method=method)

        assert md5(tempfile.name) in out_md5s, "Incorect gbk output for method {}".format(method)

