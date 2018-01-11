__version__ = "0.1.2"
import pkg_resources
try:
    version = pkg_resources.require('bioconvert')[0].version
except:
    version = __version__

import os
import shutil
import subprocess
import colorlog

# This will create a HOME/.config/bioconvert where files (e.g., executables)
# can be downloaded
from easydev import CustomConfig
configuration = CustomConfig("bioconvert", verbose=True)

if 'GOPATH' not in os.environ:
    os.environ["GOPATH"] = os.environ["HOME"]+"/go"
os.environ["PATH"] = os.environ["GOPATH"]+"/bin/:"+os.environ["PATH"]

try:
    from easydev.logging_tools import Logging
    logger = Logging("bioconvert", "INFO")
except:
    import colorlog
    logger = colorlog.getLogger("bioconvert")


def bioconvert_script(filename, where=None):
    bioconvert_path = bioconvert.__path__[0]
    share = os.path.join(bioconvert_path, 'misc')
    if where:
        filename = os.path.join(share, where, filename)
    else:
        filename = os.path.join(share, filename)
    if not os.path.exists(filename):
        raise FileNotFoundError('unknown file %s' % filename)
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
        raise FileNotFoundError('unknown file %s' % filename)
    return filename


def generate_outfile_name(infile, out_extension):
    """simple utility to replace the file extension with the given one.

    :param str infile: the path to the Input file
    :param str out_extension: Desired extension
    :return: The file path with the given extension
    :rtype: str
    """
    return '{}.{}'.format(os.path.splitext(infile)[0], out_extension)

def check_tool(executable):
    """Checks wether the given executable exists and is in the path
    
    :param str executable: name of the executable
    :return: true if the given executable exists and is in the path
    """
    return shutil.which(executable) is not None

def install_tool(executable):
    """Install the given tool, using the script:
    bioconvert/install_script/install_executable.sh
    
    :param executable to install
    :return: nothing
    """
    logger.info("Installing tool : "+executable)
    bioconvert_path = bioconvert.__path__[0]
    script = bioconvert_data('install_'+executable+'.sh', where="../scripts")
    subprocess.call(['sh',script])

import bioconvert
from bioconvert.core.base import ConvBase
from bioconvert.core.benchmark import Benchmark, BenchmarkMulticonvert
from bioconvert.core.converter import Bioconvert
from bioconvert.core.shell import shell
