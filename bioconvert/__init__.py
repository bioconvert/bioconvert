__version__ = "0.1.2"
import pkg_resources
try:
    version = pkg_resources.require('bioconvert')[0].version
except:
    version = __version__

import os
import colorlog

import bioconvert
from bioconvert.core.base import ConvBase
from bioconvert.core.benchmark import Benchmark


def init_logger():
    handler = colorlog.StreamHandler()
    formatter = colorlog.ColoredFormatter("%(log_color)s%(levelname)-8s : %(reset)s %(message)s",
                                          datefmt=None,
                                          reset=True,
                                          log_colors={
                                              'DEBUG':    'cyan',
                                              'INFO':     'green',
                                              'WARNING':  'yellow',
                                              'ERROR':    'red',
                                              'CRITICAL': 'bold_red',
                                          },
                                          secondary_log_colors={},
                                          style='%'
                                          )
    handler.setFormatter(formatter)
    logger = colorlog.getLogger('bioconvert')
    logger.addHandler(handler)
    logger.setLevel(colorlog.logging.logging.WARNING)


init_logger()


def logger_set_level(level=colorlog.logging.logging.WARNING):
    assert level in (colorlog.logging.logging.DEBUG,
                     colorlog.logging.logging.INFO,
                     colorlog.logging.logging.WARNING,
                     colorlog.logging.logging.ERROR,
                     colorlog.logging.logging.CRITICAL)
    logger = colorlog.getLogger('bioconvert')
    if level <= colorlog.logging.logging.DEBUG:
        formatter = colorlog.ColoredFormatter(
            "%(log_color)s%(levelname)-8s : %(module)s: L %(lineno)d :%(reset)s %(message)s",
            datefmt=None,
            reset=True,
            log_colors={
                'DEBUG': 'cyan',
                'INFO': 'green',
                'WARNING': 'yellow',
                'ERROR': 'red',
                'CRITICAL': 'bold_red',
            },
            secondary_log_colors={},
            style='%'
            )
        handler = logger.handlers[0]
        handler.setFormatter(formatter)

    logger.setLevel(level)


def bioconvert_data(filename, where=None):
    """
    Simple utilities to retrieve data sets from bioconvert/data directory

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
    """
    simple utility to replace the file extension with the given one.
    :param str infile: the path to the Input file
    :param str out_extension: Desired extension
    :return: The file path with the given extension
    :rtype: str
    """
    return '{}.{}'.format(os.path.splitext(infile)[0], out_extension)

