import pytest
from bioconvert.scf2fasta import Scf2Fasta, read_from_buffer, delta
from bioconvert import bioconvert_data
from easydev import TempFile, md5

def test_conv():
    # Scf V2 file
    infile_v2 = bioconvert_data("sample_v2.scf")
    expected_outfile_v2 = bioconvert_data("sample_v2.fasta")
    # Scf V3 file
    infile_v3 = bioconvert_data("sample_v3.scf")
    expected_outfile_v3 = bioconvert_data("sample_v3.fasta")

    with TempFile(suffix=".fasta") as tempfile:
        convert = Scf2Fasta(infile_v2, tempfile.name)
        convert()
        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(expected_outfile_v2)

        convert = Scf2Fasta(infile_v3, tempfile.name)
        convert()
        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(expected_outfile_v3)

def test_read_from_buffer(tmpdir):
    """Test function 'read_from_buffer(f_file, length, offset)'"""
    tmp_file = tmpdir.join("test.tmp")
    tmp_file.write(">Fake1\nWQSDESDFZQS")
    f_file = open(tmp_file)
    assert read_from_buffer(f_file, 20, 1) == "Fake1\nWQSDESDFZQS"

def test_delta():
    """Test function 'delta(rsamples, direction)'"""
    rsamples = [170, -81, -41, -25, -11, 3]
    # Bad direction
    direction = "pwet"
    with pytest.raises(Exception) as excinfo:
        delta(rsamples, direction)
    assert str(excinfo.value) == "Bad direction in 'delta'. Use\" forward\" or\" backward\"."

    # Forward
    direction = "forward"
    res = delta(rsamples, direction)
    assert res == [170, -421, 291, -24, -2, 0]

    # Backward
    direction = "backward"
    res = delta(rsamples, direction)
    assert res == [170, 259, 307, 330, 342, 357]
