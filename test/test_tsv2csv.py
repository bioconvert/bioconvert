import os

import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.tsv2csv import TSV2CSV


@pytest.mark.parametrize("method", TSV2CSV.available_methods)
def test_conv(method):
    infile = bioconvert_data("test_csv2tsv.tsv")
    expected_outile = bioconvert_data("test_csv2tsv.csv")
    with TempFile(suffix=".csv") as tempfile:
        convert = TSV2CSV(infile, tempfile.name)
        convert()
        assert md5(tempfile.name) == md5(expected_outile)
