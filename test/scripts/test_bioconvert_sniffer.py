import sys
from easydev import TempFile
import subprocess
from bioconvert import bioconvert_data
from bioconvert.scripts.sniffer import main


def test_cmdline():
    # use subprocess to prevent the image from poping up
    infile = bioconvert_data("sample_v2.scf")

    cmd = "bioconvert_sniffer  --input {}".format(infile)
    subprocess.Popen(cmd, shell=True)


def test_sniffer():
    from bioconvert.scripts.sniffer import main
    cmd = "--input {}"

    # Test SCF
    infile = bioconvert_data("sample_v2.scf")
    main(cmd.format(infile).split())


