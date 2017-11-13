from bioconvert import Benchmark
from bioconvert import bioconvert_data
from bioconvert.bam2bed import BAM2BED
from easydev import TempFile

def test_benchmark():
    input_file = bioconvert_data("test_measles.sorted.bam")
    with TempFile(suffix=".bed") as fout:
        conv = BAM2BED(input_file, fout.name)
        bench = Benchmark(conv)
        bench.plot()

        try:
            Benchmark(1)
            assert False
        except:
            assert True
