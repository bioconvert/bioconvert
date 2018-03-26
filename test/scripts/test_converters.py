import subprocess

import pytest
from easydev import TempFile, md5

from bioconvert import bioconvert_data
from bioconvert.bam2bed import BAM2BED
from bioconvert.scripts import converter


def test_converter():
    infile = bioconvert_data("test_measles.sorted.bam")
    with TempFile(suffix=".sam") as tempfile1, TempFile(suffix=".sam") as tempfile2:
        cmd = "bioconvert bam2sam %s %s --force" % (infile, tempfile1.name)
        p = subprocess.Popen(cmd, shell=True)
        assert p.wait() == 0
        import sys
        sys.argv = ["bioconvert", "bam2sam", infile, tempfile2.name, "--force"]
        converter.main()
        assert md5(tempfile1.name) == md5(tempfile2.name)


def test_converter2():
    if BAM2BED._method_bedtools.is_disabled:
        return
    infile = bioconvert_data("test_measles.sorted.bam")
    with TempFile(suffix=".bed") as tempfile:
        import sys
        sys.argv = ["bioconvert", "bam2bed", infile, tempfile.name,
                    "--method", "bedtools", "--force"]
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
                    "--output-format", "bed", "--force"]
        try:
            converter.main()
        except SystemExit:
            pass


# def test_converter_show_methods():
#     infile = bioconvert_data("test_measles.sorted.bam")
#     with TempFile(suffix=".bed") as tempfile:
#         import sys
#         sys.argv = ["bioconvert", infile, tempfile.name, "--show-methods",
#             "--force"]
#         try:
#             converter.main()
#         except SystemExit:
#             pass

def test_converter_formats():
    import sys
    sys.argv = ["bioconvert", "--formats"]
    try:
        converter.main()
        assert False
    except SystemExit:
        assert True
    except:
        assert False


def test_no_converter_specified():
    import sys
    sys.argv = ["bioconvert"]
    try:
        converter.main()
        assert False
    except SystemExit as e:
        assert e.code == 2
    except:
        assert False
    sys.argv = ["bioconvert", '--verbosity', 'DEBUG']
    with pytest.raises(Exception):
        converter.main()


def test_version():
    import sys
    sys.argv = ["bioconvert", "--version"]
    try:
        converter.main()
        assert False
    except SystemExit as e:
        assert e.code == 0
    except:
        assert False


def test_dependency_report():
    import sys
    sys.argv = ["bioconvert", "--dependency-report"]
    try:
        converter.main()
        assert False
    except SystemExit as e:
        assert e.code == 0
    except:
        assert False


def test_allow_indirect_conversion():
    import sys
    sys.argv = ["bioconvert", "-a", "--help"]
    try:
        converter.main()
        assert False
    except SystemExit as e:
        assert e.code == 0
    except:
        assert False
    sys.argv = ["bioconvert", "--allow-indirect-conversion", "--help"]
    try:
        converter.main()
        assert False
    except SystemExit as e:
        assert e.code == 0
    except:
        assert False


def test_verbose():
    infile = bioconvert_data("test_fastq2fasta_v1.fastq")
    with TempFile(suffix=".tt") as tempfile:
        import sys
        sys.argv = ["bioconvert", "-v", "CRITICAL", "fastq2fasta", infile, tempfile.name,
                    "--force"]
        converter.main()
        sys.argv = ["bioconvert", "--verbosity", "CRITICAL", "fastq2fasta", infile, tempfile.name,
                    "--force"]
        converter.main()
        sys.argv = ["bioconvert", "-l", "CRITICAL", "fastq2fasta", infile, tempfile.name,
                    "--force"]
        converter.main()
        sys.argv = ["bioconvert", "--level", "CRITICAL", "fastq2fasta", infile, tempfile.name,
                    "--force"]
        converter.main()


def test_batch():
    infile = bioconvert_data("test_fastq2fasta_v1.fastq")
    with TempFile(suffix=".tt") as tempfile:
        import sys
        sys.argv = ["bioconvert", "-m", "fastq2fasta", infile, tempfile.name,
                    "--force"]
        converter.main()


def test_converter_with_nothing():
    import sys
    sys.argv = ["bioconvert", "fastq2fasta"]
    try:
        converter.main()
        assert False
    except SystemExit as e:
        assert e.code == 2
    except:
        assert False


def test_converter_show_methods():
    import sys
    sys.argv = ["bioconvert", "fastq2fasta", "--show-methods"]
    try:
        converter.main()
    except SystemExit as e:
        assert e.code == 0


def test_indirect_conversion_without_argument():
    import sys
    infile = bioconvert_data("fastqutils_1.fastq")
    with TempFile(suffix=".clustal") as tempfile:
        sys.argv = ["bioconvert", "fastq2clustal", infile, tempfile.name, "--force"]
        # For now we want the user to explicitly indicate that (s)he agrees with an indirect conversion
        try:
            converter.main()
            assert False
        except SystemExit as e:
            assert e.code == 2
        except:
            assert False


def test_indirect_conversion():
    import sys
    infile = bioconvert_data("fastqutils_1.fastq")
    with TempFile(suffix=".clustal") as tempfile:
        sys.argv = ["bioconvert", "--allow-indirect-conversion", "-l", "DEBUG", "fastq2clustal", infile, tempfile.name,
                    "--force"]
        converter.main()
