import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.xlsx2csv import XLSX2CSV


@pytest.mark.parametrize("method", XLSX2CSV.available_methods)
def test_conv(method):
    infile = bioconvert_data("test_tabulated.xlsx")
    expected_outile = bioconvert_data("test_tabulated.csv")
    with TempFile(suffix=".csv") as tempfile:
        convert = XLSX2CSV(infile, tempfile.name)
        convert(method=method)
        assert md5(tempfile.name) == md5(expected_outile)
