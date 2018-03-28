import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.ods2csv import ODS2CSV


@pytest.mark.parametrize("method", ODS2CSV.available_methods)
def test_conv(method):
    infile = bioconvert_data("test_tabulated.ods")
    expected_outfiles = [
        bioconvert_data("test_tabulated.csv"),
        bioconvert_data("test_tabulated_with_3_more_blank_lines.csv"),
    ]
    with TempFile(suffix=".csv") as tempfile:
        convert = ODS2CSV(infile, tempfile.name)
        convert(method=method)
        assert md5(tempfile.name) in [md5(f) for f in expected_outfiles]
