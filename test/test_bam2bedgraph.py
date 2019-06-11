from bioconvert.bam2bedgraph import BAM2BEDGRAPH
from bioconvert import bioconvert_data
from easydev import TempFile, md5
import pytest

@pytest.mark.parametrize("method", BAM2BEDGRAPH.available_methods)
def test_conv(method):
    infile = bioconvert_data("test_measles.sorted.bam")
    with TempFile(suffix=".bedgraph") as tempfile:
        convert = BAM2BEDGRAPH(infile, tempfile.name)
        convert(method=method)

        # changes after 3.1 
        assert md5(tempfile.name) == "5be280e9f74e9ff1128ff1d2fe3e0812"


def test_conv_mosdepth_gz():
    infile = bioconvert_data("test_measles.sorted.bam")
    with TempFile(suffix=".bedgraph.gz") as tempfile:
        convert = BAM2BEDGRAPH(infile, tempfile.name)
        convert(method="mosdepth")


def test_conv_mosdepth_wrong_input():
    infile = bioconvert_data("test_measles.sam")


    with TempFile(suffix=".bedgraph.gz") as tempfile:
        try:
            convert = BAM2BEDGRAPH(infile, tempfile.name)
            convert(method="mosdepth")
            assert False
        except:
            False
