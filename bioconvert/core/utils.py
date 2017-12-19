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

