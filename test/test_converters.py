from easydev import TempFile
import subprocess
from bioconvert import bioconvert_data
from bioconvert.scripts import converter

def test_converter():

    infile = bioconvert_data("test_measles.sorted.bam")
    with TempFile(suffix=".bed") as tempfile:
        cmd = "bioconvert %s %s" % (infile, tempfile.name)
        subprocess.Popen(cmd, shell=True)


def test_converter1():
    infile = bioconvert_data("test_measles.sorted.bam")
    with TempFile(suffix=".bed") as tempfile:
        import sys
        sys.argv = ["bioconvert", infile, tempfile.name]
        converter.main()

