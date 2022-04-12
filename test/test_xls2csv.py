from tempfile import NamedTemporaryFile

import pytest
from easydev import TempFile, md5

from bioconvert.xls2csv import XLS2CSV

from . import test_dir


@pytest.mark.parametrize("method", XLS2CSV.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/xls/test_tabulated.xls"
    expected_outile = f"{test_dir}/data/csv/test_tabulated.csv"
    with TempFile(suffix=".csv") as tempfile:
        convert = XLS2CSV(infile, tempfile.name)
        convert(method=method)
        assert md5(tempfile.name) == md5(expected_outile)
