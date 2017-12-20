from bioconvert.bam2bigwig import BAM2BIGWIG
from bioconvert import bioconvert_data
from easydev import TempFile, md5
import pytest

def test_conv(method):
    infile = bioconvert_data('test_measles.sorted.bam')
    
    expected_outputfile = bioconvert_data('test_measles.bigwig')
    with TempFile(suffix='.bigwig') as expected_unwrapped:
        BAM2BIGWIG.unwrap_bigwig(expected_outfile, expected_unwrapped.name, strip_comment = True)
        md5out= md5(expected_unwrapped.name)

    with TempFile(suffix=".bigwig") as outfile, \
            TempFile(suffix=".bigwig") as unwrapped:
        convert = BAM2BIGWIG(infile, outfile.name)
        convert(method=method)
        BAM2BIGWIG.unwrap_fasta(
            outfile.name, unwrapped.name, strip_comment=True)
        assert md5(unwrapped.name) == md5out, \
            "{} failed".format(method)
