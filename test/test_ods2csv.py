from tempfile import NamedTemporaryFile

import pytest
from easydev import md5

from bioconvert import bioconvert_data
from bioconvert.ods2csv import ODS2CSV


@pytest.mark.parametrize("method", ODS2CSV.available_methods)
def test_conv(method):
    infile = bioconvert_data("test_tabulated.ods")
    expected_outile = bioconvert_data("test_tabulated.csv")
    with  NamedTemporaryFile(suffix=".csv", delete=False) as tempfile:
        convert = ODS2CSV(infile, tempfile.name)
        convert(method=method)
        print(tempfile.name)
        assert md5(tempfile.name) == md5(expected_outile)
