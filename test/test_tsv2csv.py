import pytest
from easydev import TempFile, md5

from bioconvert.tsv2csv import TSV2CSV

from . import test_dir

@pytest.mark.parametrize("method", TSV2CSV.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/tsv/test_tabulated.tsv"
    expected_outile = f"{test_dir}/data/csv/test_tabulated.csv"
    with TempFile(suffix=".csv") as tempfile:
        convert = TSV2CSV(infile, tempfile.name)
        convert(method=method)
        assert md5(tempfile.name) == md5(expected_outile)