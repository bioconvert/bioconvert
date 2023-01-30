import pytest
from bioconvert import TempFile, md5

from bioconvert.ods2csv import ODS2CSV

from . import test_dir


@pytest.mark.parametrize("method", ODS2CSV.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/ods/test_tabulated.ods"
    expected_outfiles = [
        f"{test_dir}/data/csv/test_tabulated.csv",
        f"{test_dir}/data/csv/test_tabulated_with_3_more_blank_lines.csv",
    ]
    with TempFile(suffix=".csv") as tempfile:
        convert = ODS2CSV(infile, tempfile.name)
        convert(method=method)
        assert md5(tempfile.name) in [md5(f) for f in expected_outfiles]
