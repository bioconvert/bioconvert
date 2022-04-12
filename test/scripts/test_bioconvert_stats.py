import sys
from easydev import TempFile
import subprocess


def test_stats():
    # use subprocess to prevent the image from poping up
    cmd = "bioconvert_stats --no-plot"
    subprocess.Popen(cmd, shell=True)


def test_stats2():
    from bioconvert.scripts.stats import main

    main(["--no-plot"])
