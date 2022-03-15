import os
import pytest
import hashlib
from easydev import TempFile, md5

from bioconvert.nexus2clustal import NEXUS2CLUSTAL

from . import test_dir

skiptravis = pytest.mark.skipif("TRAVIS_PYTHON_VERSION" in os.environ
                                and os.environ['TRAVIS_PYTHON_VERSION'].startswith("2"), 
                                reason="On travis")


@skiptravis
@pytest.mark.skipif(NEXUS2CLUSTAL._method_goalign.is_disabled, reason="missing dependencies")
def test_nexus2clustal_goalign():
    infile = f"{test_dir}/data/nexus/goalign.nexus"
    outfile = f"{test_dir}/data/clustal/goalign.clustal"
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

def test_nexus2clustal_biopython():
    infile = f"{test_dir}/data/nexus/nexus2clustal_biopython.nexus"
    outfile = f"{test_dir}/data/clustal/nexus2clustal_biopython.clustal"
    with TempFile(suffix=".nexus") as tempfile:
        converter = NEXUS2CLUSTAL(infile, tempfile.name)
        converter(method='biopython')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)

def test_nexus2clustal_squizz():
    infile = f"{test_dir}/data/nexus/nexus2clustal_squizz.nexus"
    outfile = f"{test_dir}/data/clustal/nexus2clustal_squizz.clustal"
    with TempFile(suffix=".nexus") as tempfile:
        converter = NEXUS2CLUSTAL(infile, tempfile.name)
        converter(method='squizz')

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)