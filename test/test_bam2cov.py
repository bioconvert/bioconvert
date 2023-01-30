from bioconvert.bam2cov import BAM2COV
from bioconvert import TempFile, md5
import pytest

from . import test_dir


@pytest.mark.parametrize("method", BAM2COV.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/bam/test_measles.sorted.bam"
    with TempFile(suffix=".cov") as tempfile:
        convert = BAM2COV(infile, tempfile.name)
        convert(method=method)

        # Check that the output is correct with a checksum
        # Note that we cannot test the md5 on a gzip file but only
        # on the original data. This check sum was computed
        # fro the unzipped version of biokit/data/converters/measles.bed
        assert md5(tempfile.name) == "84702e19ba3a27900f271990e0cc72a0"
