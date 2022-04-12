import gzip
import bz2

import pytest

from bioconvert.gz2bz2 import GZ2BZ2
from easydev import TempFile


@pytest.mark.parametrize("method", GZ2BZ2.available_methods)
def test_conv(method):
    # generate a fake sequence
    content_ref = b"atgc" * 50
    # Compress it
    gzip_string = gzip.compress(content_ref, 9)
    # Generate two tmp file
    gzip_file = TempFile(suffix="gzip")
    gz_file = TempFile(suffix="gz")
    # Write the gzip file with compressed fake sequence
    with open(gzip_file.name, "wb") as gzip_to_write:
        gzip_to_write.write(gzip_string)

    # Convert gzip to gz
    converter = GZ2BZ2(gzip_file.name, gz_file.name)
    converter(method=method)

    # gunzip result file
    with bz2.open(gz_file.name, "rb") as gz_to_read:
        content = gz_to_read.read()

    # check conversion
    assert content == content_ref
