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

        # Check that the output is correct with a checksum
        # Note that we cannot test the md5 on a gzip file but only 
        # on the original data. This check sum was computed
        # from the bedgraph file computed with bedtools
        # from  test_measles.sorted.bam 
        assert md5(tempfile.name) in ["86d7a2236c40ddb4c3a8dccdbd978356", "c0e066803991efdd25af5fa7d175a9df"]
        #"c0e066803991efdd25af5fa7d175a9df" is for ucsc.bedgraph but this is not the integer format 
        # returned by bedtools for now because no normalization is applied

        #convert = BAM2BEDGRAPH(infile, tempfile.name)
        #convert(method="bedtools")
        #assert md5(tempfile.name) == "c0e066803991efdd25af5fa7d175a9df"
