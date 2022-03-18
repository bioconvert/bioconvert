import subprocess
import os
import shutil
import pytest
from easydev import TempFile, md5
import sys

from bioconvert.bam2cov import BAM2COV
from bioconvert.scripts import converter

from . import test_dir

fastq_file = f"{test_dir}/../data/fastq/test_fastq2fasta_v1.fastq"
fasta_file = f"{test_dir}/../data/fasta/test_fastq2fasta_v1.fasta"

def test_converter_compression():

    with TempFile(suffix=".fasta.gz") as tempfile:
        sys.argv = ["bioconvert", "fastq2fasta", fastq_file, tempfile.name, "--force"]
        converter.main()

    with TempFile(suffix=".fasta.bz2") as tempfile:
        sys.argv = ["bioconvert", "fastq2fasta", fastq_file, tempfile.name, "--force"]
        converter.main()

    with TempFile(suffix=".fastq.gz") as tempfile:
        sys.argv = ["bioconvert", "fasta2fastq", fasta_file, tempfile.name, "--force"]
        converter.main()

    with TempFile(suffix=".fastq.bz2") as tempfile:
        sys.argv = ["bioconvert", "fasta2fastq", fasta_file, tempfile.name, "--force"]
        converter.main()

    with TempFile(suffix=".fastq.dsrc") as tempfile:
        sys.argv = ["bioconvert", "fasta2fastq", fasta_file, tempfile.name, "--force"]
        converter.main()


def test_converter_wrong_input_file():
    with TempFile(suffix=".fasta") as tempfile1, TempFile(suffix=".fasta") as tempfile2:
        cmd = "bioconvert fastq2fasta {} {} --force".format("missing.fastq", tempfile1.name)
        p = subprocess.Popen(cmd, shell=True)
        assert p.wait() == 1
        sys.argv = ["bioconvert", "fastq2fasta", "missing.fastq", tempfile2.name, "--force"]
        try:
            converter.main()
            assert False
        except:
            assert True


def test_converter():
    infile = fastq_file
    with TempFile(suffix=".fasta") as tempfile1, TempFile(suffix=".fasta") as tempfile2:
        cmd = "bioconvert fastq2fasta {} {} --force".format(infile, tempfile1.name)
        p = subprocess.Popen(cmd, shell=True)
        assert p.wait() == 0
        sys.argv = ["bioconvert", "fastq2fasta", infile, tempfile2.name, "--force"]
        converter.main()
        assert md5(tempfile1.name) == md5(tempfile2.name)


def test_converter_without_converter():
    infile = fastq_file
    with TempFile(suffix=".fasta") as tempfile1, TempFile(suffix=".fasta") as tempfile2:
        cmd = "bioconvert {} {} --force".format(infile, tempfile1.name)
        p = subprocess.Popen(cmd, shell=True)
        assert p.wait() == 0
        sys.argv = ["bioconvert", infile, tempfile2.name, "--force"]
        converter.main()
        assert md5(tempfile1.name) == md5(tempfile2.name)


@pytest.mark.skipif(BAM2COV._method_bedtools.is_disabled, reason="missing dependencies")
def test_converter2():
    infile = f"{test_dir}/../data/bam/test_measles.sorted.bam"
    with TempFile(suffix=".bed") as tempfile:
        sys.argv = ["bioconvert", "bam2cov", infile, tempfile.name,
                    "--method", "bedtools", "--force"]
        converter.main()


# FIXME 3/3/22
def _test_plink_no_extension():

    infile = f"{test_dir}/../data/plink/plink_toy.ped"
    infile = infile.replace(".ped", "")

    with TempFile(suffix="") as outfile:
        sys.argv = ["bioconvert", "plink2bplink", infile, outfile.name, "--force"]
        converter.main()


def test_converter_no_outfile():
    infile = fastq_file
    with TempFile(suffix=".fastq") as tempfile:
        shutil.copy(infile, tempfile.name)
        sys.argv = ["bioconvert", "fastq2fasta", tempfile.name, "--force", "--raise-exception"]
        converter.main()
        # os.remove(infile[:-3] + "sam")


def test_converter_no_outfile_without_srs_argv():
    infile = fastq_file
    with TempFile(suffix=".fastq") as tempfile:
        shutil.copy(infile, tempfile.name)
        args = ["bioconvert", "fastq2fasta", tempfile.name, "--force", "--raise-exception"]
        converter.main(args[1:])


def test_converter_output_format():
    infile = f"{test_dir}/../data/bam/test_measles.sorted.bam"
    with TempFile() as tempfile:
        sys.argv = ["bioconvert", infile, tempfile.name,
                    "--output-format", "bed", "--force"]
        try:
            converter.main()
        except SystemExit:
            pass


def test_no_converter_specified():
    try:
        sys.argv = ["bioconvert"]
        converter.main()
        assert False
    except SystemExit as e:
        assert e.code == 1
    except Exception as e:
        assert type(e) == SystemExit
    try:
        sys.argv = ["bioconvert", '--verbosity', 'DEBUG']
        converter.main()
        assert False
    except SystemExit as e:
        assert e.code == 2
    except Exception as e:
        assert type(e) == SystemExit and e.code == 2


def test_not_existing_param():
    sys.argv = ["bioconvert", "--tagada"]
    try:
        converter.main()
        assert False
    except SystemExit as e:
        assert e.code == 2
    except:
        assert False


def test_not_existing_subcommand():
    sys.argv = ["bioconvert", "bam2tagada", "--help"]
    try:
        converter.main()
        assert False
    except SystemExit as e:
        assert e.code in [1,2]
    except:
        assert False


def test_version():
    sys.argv = ["bioconvert", "--version"]
    try:
        converter.main()
        assert False
    except SystemExit as e:
        assert e.code == 0
    except:
        assert False


def test_help_works():
    sys.argv = ["bioconvert", "--help"]
    try:
        converter.main()
        assert False
    except SystemExit as e:
        assert e.code == 0
    except:
        assert False
    sys.argv = ["bioconvert", "-v", "DEBUG", "--help"]
    try:
        converter.main()
        assert False
    except SystemExit as e:
        assert e.code == 0
    except:
        assert False
    sys.argv = ["bioconvert", "--help", "-v", "DEBUG"]
    try:
        converter.main()
        assert False
    except SystemExit as e:
        assert e.code == 0
    except:
        assert False
    sys.argv = ["bioconvert", "fastq2fasta", "--help"]
    try:
        converter.main()
        assert False
    except SystemExit as e:
        assert e.code == 0
    except:
        assert False
    sys.argv = ["bioconvert", "fastq2fasta", "--help", "-v", "DEBUG"]
    try:
        converter.main()
        assert False
    except SystemExit as e:
        assert e.code == 0
    except:
        assert False


def test_dependency_report():
    sys.argv = ["bioconvert", "--dependency-report"]
    try:
        converter.main()
        assert False
    except SystemExit as e:
        assert e.code == 0
    except:
        assert False


def test_allow_indirect_conversion():
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


def test_wrong_verbose():
    sys.argv = ["bioconvert", "--version", "-v" "INFO "]
    try:
        converter.main()
        assert False
    except SystemExit as e:
        assert e.code == 2
    except:
        assert False
    sys.argv = ["bioconvert", "--version", "-v" "TRALALA"]
    try:
        converter.main()
        assert False
    except SystemExit as e:
        assert e.code == 2
    except:
        assert False
    sys.argv = ["bioconvert", "fastq2fasta", "--show-methods", "-v", "TRALALA"]
    try:
        converter.main()
        assert False
    except SystemExit as e:
        assert e.code == 2
    except:
        assert False


def test_verbose():
    infile = fastq_file
    with TempFile(suffix=".tt") as tempfile:
        sys.argv = ["bioconvert", "fastq2fasta", infile, tempfile.name,
                    "--force", "-v", "CRITICAL"]
        converter.main()
        sys.argv = ["bioconvert","fastq2fasta", infile, tempfile.name,
                    "--force","--verbosity", "CRITICAL"]
        converter.main()


def test_close_match():
    infile = fastq_file
    with TempFile(suffix=".tt") as tempfile:
        sys.argv = ["bioconvert", "fastq2fastpp", infile, tempfile.name,
                    "--force"]
        try:
            converter.main()
            assert False
        except SystemExit as e:
            assert e.code == 1
        except:
            assert False


def test_batch():
    infile = fastq_file
    with TempFile(suffix=".tt") as tempfile:
        sys.argv = ["bioconvert", "fastq2fasta", infile, tempfile.name,
                    "--force", "-m"]
        converter.main()


def test_converter_with_nothing():
    sys.argv = ["bioconvert", "fastq2fasta"]
    try:
        converter.main()
        assert False
    except SystemExit as e:
        assert e.code >0
    except:
        assert False


def test_converter_show_methods():
    sys.argv = ["bioconvert", "fastq2fasta", "--show-methods"]
    try:
        converter.main()
    except SystemExit as e:
        assert e.code == 0

def test_indirect_conversion_without_argument():
    infile = f"{test_dir}/../data/fastq/ERR3295124.fastq"
    with TempFile(suffix=".clustal") as tempfile:
        sys.argv = ["bioconvert", "fastq2clustal", infile, tempfile.name, "--force"]
        # For now we want the user to explicitly indicate that (s)he agrees 
        # with an indirect conversion
        try:
            converter.main()
            assert False
        except SystemExit as e:
            assert e.code >0
        except:
            assert False


def test_indirect_conversion():
    infile = f"{test_dir}/../data/fastq/ERR3295124.fastq"
    #infile = bioconvert_data("fastqutils_1.fastq")
    with TempFile(suffix=".clustal") as tempfile:
        sys.argv = ["bioconvert", "fastq2clustal", infile, tempfile.name,
                    "--force", "--allow-indirect-conversion", "-v", "DEBUG"]
        converter.main()
        sys.argv = ["bioconvert", "fastq2clustal", infile, tempfile.name,
                    "--force", "-a", "-v", "DEBUG"]
        converter.main()


def test_conversion_graph_error():
    sys.argv = ["bioconvert", "--conversion-graph", "toto"]
    try:
        converter.main()
    except SystemExit as e:
        assert e.code >0


def test_conversion_graph():
    sys.argv = ["bioconvert", "--conversion-graph", "cytoscape"]
    try:
        converter.main()
    except SystemExit as e:
        assert e.code == 0
    sys.argv = ["bioconvert", "--conversion-graph", "cytoscape-all"]
    try:
        converter.main()
    except SystemExit as e:
        assert e.code == 0


def is_osx():
    if "TRAVIS_OS_NAME" in os.environ:
        if os.environ["TRAVIS_OS_NAME"] == "osx":
            return True
    return False


@pytest.mark.skipif("DISPLAY" not in os.environ, reason="no DISPLAY available, will fail otherwise")
@pytest.mark.skipif(is_osx(), reason="unknown failure on travis april 2019")
def test_converter_benchmark():
    infile = fastq_file
    with TempFile(suffix=".fasta") as tempfile1, TempFile(suffix=".fasta") as tempfile2:
        cmd = "bioconvert fastq2fasta {} {} --force -b".format(infile, tempfile1.name)
        p = subprocess.Popen(cmd, shell=True)
        assert p.wait() == 0
        sys.argv = ["bioconvert", "fastq2fasta", infile, tempfile2.name, "--force", "-b"]
        converter.main()
        assert md5(tempfile1.name) == md5(tempfile2.name)