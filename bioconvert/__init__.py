__version__ = "0.1.2"
import pkg_resources
try:
    version = pkg_resources.require(bioconvert)[0].version
except:
    version = __version__

# Creates the data directory if it does not exist
from easydev import CustomConfig
PATH = CustomConfig("bioconvert").user_config_dir

import colorlog as logger
def bioconvert_debug_level(level="WARNING"):
    """A deubg level setter at top level of the library"""
    assert level in ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
    logging_level = getattr(logger.logging.logging, level)
    logger.getLogger().setLevel(logging_level)

def bioconvert_data(filename, where=None):
    """Simple utilities to retrieve data sets from bioconvert/data directory"""
    import os
    import easydev
    bioconvert_path = easydev.get_package_location('bioconvert')
    share = os.sep.join([bioconvert_path , "bioconvert", 'data'])
    # in the code one may use / or \ 
    if where:
        filename = os.sep.join([share, where, filename])
    else:
        filename = os.sep.join([share, filename])
    if os.path.exists(filename) is False:
        raise Exception('unknown file %s' % filename)
    return filename


from bioconvert.core.base import ConvBase
