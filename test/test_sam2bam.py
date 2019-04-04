import pytest
from easydev import TempFile, md5
import filecmp 

from bioconvert import bioconvert_data
from bioconvert.sam2bam import SAM2BAM
from bioconvert.bam2sam import BAM2SAM


@pytest.mark.skipif(len(SAM2BAM.available_methods) == 0, reason="missing dependencies")
def test_conv():
    infile = bioconvert_data("test_measles.sam")
    outfile = bioconvert_data("test_measles.bam")
    with TempFile(suffix=".bam") as tempfile:
        convert = SAM2BAM(infile, tempfile.name)
        convert()

        # Difficult to test the md5 since bam is binary 
        # recurrent problem is the version stored in the file that keeps
        # changing
        # However, we can reverse back the bam 2 sam and they should agree since
        # the same does not store the version
        with TempFile(suffix=".sam") as samfile:
            convert = BAM2SAM(tempfile.name, samfile.name)
            convert()
            filecmp.cmp(samfile.name, infile)