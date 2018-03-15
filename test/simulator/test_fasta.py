from bioconvert.simulator import fasta
from easydev import TempFile


def test_fasta():
    with TempFile(suffix=".fasta") as fout:
        f = fasta.FastaSim(fout.name)
        f.nreads = 1000
        f.simulate()
