from bioconvert.simulator import fastq
from bioconvert import TempFile


def test_fastq():

    with TempFile(suffix=".fastq") as fout:
        f = fastq.FastqSim(fout.name)
        f.nreads = 1000
        f.simulate()
