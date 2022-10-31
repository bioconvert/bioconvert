import pytest
from easydev import TempFile, md5

from bioconvert.gff32gtf import GFF32GTF

from . import test_dir


@pytest.mark.parametrize("method", GFF32GTF.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/GFF3/gff3_example.gff"
    with TempFile(suffix=".tsv") as tempfile:
        convert = GFF32GTF(infile, tempfile.name)
        convert(method=method)
