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


def test_converter2():
    infile = bioconvert_data("test_measles.sorted.bam")
    with TempFile(suffix=".bed") as tempfile:
        import sys
        sys.argv = ["bioconvert", infile, tempfile.name, 
                    "--method" , "bedtools"]
        converter.main()


def test_converter_no_outfile():
    infile = bioconvert_data("test_measles.sorted.bam")
    try:
        sys.argv = ["bioconvert", infile]
        converter.main()
    except:
        pass


def test_converter_no_infile_ext():
    try:
        sys.argv = ["bioconvert", "test_without_ext", "--input-format", "bam"]
        converter.main()
    except:
        pass

def test_converter_output_format():
    infile = bioconvert_data("test_measles.sorted.bam")
    with TempFile() as tempfile:
        import sys
        sys.argv = ["bioconvert", infile, tempfile.name,
                    "--output-format" , "bed"]
        try:
            converter.main()
        except SystemExit:
            pass

def test_converter_show_methods():
    infile = bioconvert_data("test_measles.sorted.bam")
    with TempFile(suffix=".bed") as tempfile:
        import sys
        sys.argv = ["bioconvert", infile, tempfile.name, "--show-methods"]
        try:
            converter.main()
        except SystemExit:
            pass

def test_converter_formats():
    infile = bioconvert_data("test_measles.sorted.bam")
    with TempFile(suffix=".bed") as tempfile:
        import sys
        sys.argv = ["bioconvert", "--formats"]
        try:
            converter.main()
            assert False
        except SystemExit:
            assert True
        except:
            assert False

