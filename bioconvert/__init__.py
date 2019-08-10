__version__ = "0.1.2"
import pkg_resources
try:
    version = pkg_resources.require('bioconvert')[0].version
except:
    version = __version__

import os
import colorlog

# This will create a HOME/.config/bioconvert where files (e.g., executables)
# can be downloaded
from easydev import CustomConfig
configuration = CustomConfig("bioconvert", verbose=True)

os.environ["GOPATH"]= os.environ["HOME"]+"/go"
os.environ["PATH"] = os.environ["GOPATH"]+"/bin/:"+os.environ["PATH"]

from easydev.logging_tools import Logging
logger = Logging("bioconvert", "INFO")


def bioconvert_script(filename, where=None):
    bioconvert_path = bioconvert.__path__[0]
    share = os.path.join(bioconvert_path, 'misc')
    if where:
        filename = os.path.join(share, where, filename)
    else:
        filename = os.path.join(share, filename)
    if not os.path.exists(filename):
        raise FileNotFoundError('unknown file {}'.format(filename))
    return filename


def bioconvert_data(filename, where=None):
    """Simple utilities to retrieve data sets from bioconvert/data directory

    :param str filename: the name of the data file to get the path
    :param str where:
    """
    bioconvert_path = bioconvert.__path__[0]
    share = os.path.join(bioconvert_path, 'data')
    # in the code one may use / or \ 
    if where:
        filename = os.path.join(share, where, filename)
    else:
        filename = os.path.join(share, filename)
    if not os.path.exists(filename):
        raise FileNotFoundError('unknown file {}'.format(filename))
    return filename


def info():
    from bioconvert.core.registry import Registry
    r = Registry()
    info = r.get_info()
    converters = [x for x in info.items()]
    data = [info[k] for k,v in info.items()]
    msg = "Bioconvert contains {} converters including {} methods"
    return msg.format(len(converters), sum(data)) 

import bioconvert
from bioconvert.core.base import ConvBase
from bioconvert.core.registry import Registry
from bioconvert.core.decorators import requires, compressor
from bioconvert.core.benchmark import Benchmark, BenchmarkMulticonvert
from bioconvert.core.converter import Bioconvert
from bioconvert.core.shell import shell
from bioconvert.io.sniffer import Sniffer
