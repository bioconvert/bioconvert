# -*- coding: utf-8 -*-
#
#  This file is part of Bioconvert software
#
#  Copyright (c) 2017 - Bioconvert Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/biokit/bioconvert
#  documentation: http://bioconvert.readthedocs.io
#
##############################################################################
"""misc utility functions """
import os
import sys
import bioconvert
from bioconvert.core import extensions

def get_extension(filename, remove_compression=False):
    """Return extension of a filename

    ::

        >>> get_extension("test.fastq")
        fastq
        >>> get_extension("test.fastq.gz")
        fastq

    """
    compression = [".gz", ".bz2", ".bzip2", ".dsrc"]
    if remove_compression is True:
        # remove the .gz, .bz2, ..., extensions
        for this in compression:
            if filename.endswith(this):
                filename = filename.rsplit(".", 1)[0]

    if len(filename.split(".")) == 1:
        return None
    else:
        return os.path.splitext(filename)[-1]


def generate_outfile_name(infile, out_extension):
    """simple utility to replace the file extension with the given one.

    :param str infile: the path to the Input file
    :param str out_extension: Desired extension
    :return: The file path with the given extension
    :rtype: str
    """
    return '{}.{}'.format(os.path.splitext(infile)[0], out_extension)

def get_format_from_extension(extension):
    """get format from extension.

    :param extension: the extension
    :return: the corresponding format
    :rtype: str
    """
    extension = extension.lstrip(".")
    for fmt, ext_list in extensions.items():
        if extension in ext_list:
            return fmt.upper()

    # The extension was not found
    bioconvert.logger.critical("No format was found for extension '%s'" % extension)
    bioconvert.logger.critical("Use --formats to know the available formats and --help"
                               " for examples")
    sys.exit(1)

