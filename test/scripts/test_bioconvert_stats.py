import sys
from easydev import TempFile
import subprocess
from bioconvert.scripts.init_convert import main

def test_init_converter():
    cmd = "bioconvert_stats --no-plot"
    subprocess.Popen(cmd, shell=True)

#how to prevent the image from poping up
#def test_init_converter():
#    cmd = "bioconvert_stats"
#    subprocess.Popen(cmd, shell=False)



