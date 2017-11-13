from bioconvert import Bioconvert
from bioconvert import bioconvert_data
from bioconvert.bam2bed import BAM2BED
from easydev import TempFile

def test_bioconvert():
    infile = bioconvert_data("test_measles.sorted.bam")
    with TempFile(suffix=".fasta") as fout:
        c = Bioconvert(infile, fout.name)
        c()
        c.boxplot_benchmark()



