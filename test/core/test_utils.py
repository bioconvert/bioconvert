from bioconvert.core import utils


83 - 86, 93 - 94, 96 - 97, 99 - 100


def test_utils():

    assert utils.get_extension("test.fastq") == "fastq"
    assert utils.get_extension("test.fastq.gz", remove_compression=True) == "fastq"
    assert utils.get_extension("test") is None

    assert utils.generate_outfile_name("test.fastq", "fq") == "test.fq"

    assert utils.get_format_from_extension(".bam") == "BAM"
    assert utils.get_format_from_extension(".bb") == "BIGBED"
    try:
        utils.get_format_from_extension(".temp")
        assert False
    except:
        assert True
