from bioconvert.scf2fastq import Scf2Fastq
from bioconvert import bioconvert_data
from easydev import TempFile, md5

def test_conv():
    # Scf V2 file
    infile_v2 = bioconvert_data("sample_v2.scf")
    expected_outfile_v2 = bioconvert_data("sample_v2.fastq")
    # Scf V3 file
    infile_v3 = bioconvert_data("sample_v3.scf")
    expected_outfile_v3 = bioconvert_data("sample_v3.fastq")

    with TempFile(suffix=".fastq") as tempfile:
        convert = Scf2Fastq(infile_v2, tempfile.name)
        convert()
        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(expected_outfile_v2)

        convert = Scf2Fastq(infile_v3, tempfile.name)
        convert()
        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(expected_outfile_v3)
