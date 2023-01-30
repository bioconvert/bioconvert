import pytest
from bioconvert import TempFile, md5

from bioconvert.xlsx2csv import XLSX2CSV

from . import test_dir


@pytest.mark.parametrize("method", XLSX2CSV.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/xlsx/test_tabulated.xlsx"
    expected_outile = f"{test_dir}/data/csv/test_tabulated.csv"
    with TempFile(suffix=".csv") as tempfile:
        convert = XLSX2CSV(infile, tempfile.name)
        convert(method=method)
        assert md5(tempfile.name) == md5(expected_outile)
