
from bioconvert.core import utils



def test_utils():

    assert utils.get_extension("test.fastq") == ".fastq"
