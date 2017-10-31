from bioconvert import Benchmark
from bioconvert import bioconvert_data
from bioconvert.bam2bed import BAM2BED


def test_benchmark():
    input_file = bioconvert_data("test_measles.sorted.bam")
    conv = BAM2BED(input_file, "test.bed")
    bench = Benchmark(conv)
    bench.plot()

    try:
        Benchmark(1)
        assert False
    except:
        assert True
