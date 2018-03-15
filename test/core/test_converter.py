from bioconvert import Bioconvert
from bioconvert import bioconvert_data
from easydev import TempFile


def test_bioconvert():
    infile = bioconvert_data("test_measles.sorted.bam")
    with TempFile(suffix=".fasta") as fout:
        c = Bioconvert(infile, fout.name, force=True)
        c()
        c.boxplot_benchmark()


def test_indirect_conversion():
    infile = bioconvert_data("fastqutils.sam")
    with TempFile(suffix=".fasta") as fout:
        c = Bioconvert(infile, fout.name, force=True)
        c()
        c.boxplot_benchmark()
