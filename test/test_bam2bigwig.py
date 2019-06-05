from bioconvert.bam2bigwig import BAM2BIGWIG
from bioconvert.bam2bedgraph import BAM2BEDGRAPH
from bioconvert import bioconvert_data
from easydev import TempFile, md5
import pytest
import os

# commented due to constant failure on travis with py3.5


# Here we will scan all available methods and repeat the test
# automatically for each method


skiptravis = pytest.mark.skipif( "TRAVIS_PYTHON_VERSION" in os.environ,
    reason="fails on travis (deeptools and numpy not compatible)")



@skiptravis
@pytest.mark.parametrize("method", BAM2BIGWIG.available_methods)
def _test_conv(method):

    # the input file
    infile = bioconvert_data('test_measles.sorted.bam')

    # What is the expected md5sum of the final output file ?
    expected_outputfile = bioconvert_data('test_measles.bigwig')
    md5out= md5(expected_outputfile)

    # Call convert and check that the output file created has the correct md5sum
    with TempFile(suffix=".bigwig") as outfile:
        convert = BAM2BIGWIG(infile, outfile.name)
        if (method == 'ucsc'):
            convert(method=method, chrom_sizes=bioconvert_data("hg38.chrom.sizes"))
            # TODO Failed in OCt 2018 . why ? bamCoverage vesrion in header ?
            #assert md5(outfile.name) == '61abd0de51bd614136ad85ae0a1ff85b', "{} failed".format(method)
        else:
            convert(method=method)
            # TODO. Failed in oct 2018. why . bamCoverage version in header ?
            #assert md5(outfile.name) == md5out, "{} failed".format(method)
