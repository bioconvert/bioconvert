import os

import pytest
from easydev import TempFile, md5

from bioconvert.csv2tsv import CSV2TSV

from . import test_dir


@pytest.mark.parametrize("method", CSV2TSV.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/csv/test_tabulated.csv"
    expected_outile = f"{test_dir}/data/tsv/test_tabulated.tsv"
    with TempFile(suffix=".tsv") as tempfile:
        convert = CSV2TSV(infile, tempfile.name)
        convert(method=method)
        assert md5(tempfile.name) == md5(expected_outile)
