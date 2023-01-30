from tempfile import NamedTemporaryFile

import pytest
from bioconvert import TempFile, md5

from bioconvert.csv2xls import CSV2XLS
from bioconvert.xls2csv import XLS2CSV

from . import test_dir


@pytest.mark.parametrize("method", CSV2XLS.available_methods)
def test_conv(method):
    # XLS file may contains bold, border, ... i then prefere to convert
    # it back to csv to check if it is ok or not

    infile = f"{test_dir}/data/csv/test_tabulated.csv"
    expected_outile = f"{test_dir}/data/csv/test_tabulated.csv"

    with TempFile(suffix=".csv") as temp_csv, TempFile(suffix=".xls") as temp_xls:
        convert = CSV2XLS(infile, temp_xls.name)
        convert(method=method)
        convert = XLS2CSV(temp_xls.name, temp_csv.name)
        convert(method=method)
        assert md5(temp_csv.name) == md5(expected_outile)
