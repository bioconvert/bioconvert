import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.nexus2clustal import NEXUS2CLUSTAL

import re


@pytest.mark.parametrize("method", NEXUS2CLUSTAL.available_methods)
def test_nx2aln(method):
    infile = bioconvert_data(method + ".nexus")
    outfile = bioconvert_data(method + ".clustal")
    with TempFile(suffix=".clustal") as tempfile:
        converter = NEXUS2CLUSTAL(infile, tempfile.name)
        converter(method=method)
        with TempFile(suffix=".clustal") as clean:
            with open(clean.name, 'w') as cleant:
                with open(tempfile.name) as infile:
                    for line in infile:
                        ll = re.sub(
                            r"^CLUSTAL.*",
                            r"CLUSTAL",
                            line
                        )
                        cleant.write(ll)
                        print(ll)
            # Check that the output is correct with a checksum
            assert md5(clean.name) == md5(outfile)
