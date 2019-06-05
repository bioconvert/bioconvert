from bioconvert import Sniffer, bioconvert_data
import os
import glob
import bioconvert
import pytest

# let us test all formats
from bioconvert.core.extensions import extensions
formats = list(extensions.keys())

bioconvert_path = bioconvert.__path__[0]
share = os.path.join(bioconvert_path, 'data')

#os.environ['PYTHONHASHSEED'] = "0"

s = Sniffer()

def test_sniffer():
    # Sniff file not included
    filename = glob.glob(share + "/__init__.py")[0]
    assert s.sniff(filename) is None

    # now with a valid file already implemented
    assert s.sniff(bioconvert_data("biopython.clustal")) == "clustal"


# here we loop cross al data files provided in bioconvert.
# Ideally, the sniffer should include all formats.
@pytest.mark.parametrize("frmt", formats)
def test_sniffer_all(frmt):
    exts = extensions[frmt]
    for ext in exts:
        filenames = glob.glob(share + "/*{}".format(ext))
        for filename in filenames:
            ret = s.sniff(filename)

            # ambiguity with files ending in gff that can be either gff2 or gff3
            # we can 
            #if ext == "gff":
            #    assert ret in ["gff2", "gff3"]
            #else:
            #    assert ret in [frmt, None], "{}; ret:{} frmt:{}".format(ret, frmt, filename)


