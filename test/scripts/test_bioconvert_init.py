import sys
from easydev import TempFile
import subprocess
from bioconvert.scripts.init_convert import main


def test_init_converter():
    cmd = "bioconvert_init -i bam -o sam"
    subprocess.Popen(cmd, shell=True)


def test_init_converter2():
    sys.argv = ["bioconvert_init", "-i", "bam", "-o", "sam"]
    main()
