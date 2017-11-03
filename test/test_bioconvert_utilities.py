import os
import pytest
import bioconvert
from bioconvert import generate_outfile_name, bioconvert_data


def test_bioconvert_data():
    file_name = 'squizz.phylip'
    assert bioconvert_data(file_name) == os.path.join(bioconvert.__path__[0], 'data', file_name)
    with pytest.raises(FileNotFoundError, message="Excepting FileNotFoundError"):
        file_name = 'foo.bar'
        assert bioconvert_data(file_name) == os.path.join(bioconvert.__path__[0], 'data', file_name)




def test_generate_outfile_name():
    assert generate_outfile_name('foo.fasta', 'phylip') == 'foo.phylip'
    assert generate_outfile_name('/foo/bar.ext.fasta', 'clustal') == '/foo/bar.ext.clustal'


