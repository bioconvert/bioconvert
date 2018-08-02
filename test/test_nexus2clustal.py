import os
import pytest
import hashlib
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.nexus2clustal import NEXUS2CLUSTAL

skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ
                                and os.environ['TRAVIS_PYTHON_VERSION'].startswith("2"), 
                                reason="On travis")


@skiptravis
def test_nexus2clustal_goalign():
    infile = bioconvert_data("goalign.nexus")
    outfile = bioconvert_data("goalign.clustal")
    with TempFile(suffix=".nexus") as tempfile:
        converter = NEXUS2CLUSTAL(infile, tempfile.name)
        converter(method='goalign')

        ## We remove goalign version from the first line
        out = ""
        with open(tempfile.name) as f:
            lines = f.readlines()
            if len(lines)>0:
                clustal = lines[0].split(" ")
                if len(clustal) > 0:
                    lines[0]=clustal[0]+"\n"
            out = ''.join(lines)

        # Check that the output is correct with a checksum
        assert hashlib.md5(out.encode('utf-8')).hexdigest() == md5(outfile)
