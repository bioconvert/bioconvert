import os
import pytest
import bioconvert
from bioconvert import generate_outfile_name, bioconvert_data


def test_bioconvert_data():
    import sys
    file_name = 'squizz.phylip'
    print("####### interpreter", sys.executable, file=sys.stderr)
    print("####### bioconvert_data(file_name)", bioconvert_data(file_name), file=sys.stderr)
    print("####### bioconvert.__path__[0]    ", os.path.join(bioconvert.__path__[0], 'data', file_name), file=sys.stderr)
    print("####### sys.path                  \n",
          '\n'.join(['#######                           {}'.format(p) for p in sys.path]), sep='')
    sys.stderr.flush()
    raise Exception()
    assert bioconvert_data(file_name) == os.path.join(bioconvert.__path__[0], 'data', file_name)
    with pytest.raises(FileNotFoundError, message="Excepting FileNotFoundError"):
        file_name = 'foo.bar'
        assert bioconvert_data(file_name) == os.path.join(bioconvert.__path__[0], 'data', file_name)




def test_generate_outfile_name():
    assert generate_outfile_name('foo.fasta', 'phylip') == 'foo.phylip'
    assert generate_outfile_name('/foo/bar.ext.fasta', 'clustal') == '/foo/bar.ext.clustal'


