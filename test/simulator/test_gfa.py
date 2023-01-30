from bioconvert.simulator import gfa
from bioconvert import TempFile


def test_gfa():
    with TempFile(suffix=".gfa") as fout:
        f = gfa.GFASim(fout.name)
        f.simulate()
