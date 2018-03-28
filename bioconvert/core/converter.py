# -*- coding: utf-8 -*-

##############################################################################
#  This file is part of Biokit software
#
#  Copyright (c) 2017 - Biokit Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/biokit/bioconvert
#  documentation: http://bioconvert.readthedocs.io
##############################################################################
""".. rubric:: Standalone application dedicated to conversion"""
import os

import colorlog

from bioconvert.core.base import ConvMeta

_log = colorlog.getLogger(__name__)

from bioconvert.core.base import make_chain
from bioconvert.core.registry import Registry
from bioconvert.core.utils import get_extension as getext
from bioconvert.core.utils import get_format_from_extension


__all__ = ['Bioconvert']


class Bioconvert(object):
    """Universal converter used by the standalone

    ::

        from bioconvert import Bioconvert
        c = Bioconvert("test.fastq", "test.fasta")


    """
    def __init__(self, infile, outfile, in_fmt=None, out_fmt=None, force=False):
        """.. rubric:: constructor

        :param str infile: The path of the input file.
        :param str outfile: The path of The output file
        :param str in_fmt: the format for the input file
        :param str out_fmt: the format for the output
        :param bool force: overwrite output file if it exists already
            otherwise raises an error

        """
        # don't check the input file because there are cases where input parameter is just a prefix
        # if os.path.exists(infile) is False:
        #     msg = "Incorrect input file: %s" % infile
        #     _log.error(msg)
        #     raise ValueError(msg)

        # check existence of output file. If it exists,
        # fails except if force argument is set to True
        if os.path.exists(outfile) is True:
            msg = "output file {} exists already".format(outfile)
            _log.warning("output file exists already")
            if force is False:
                _log.critical("output file exists. If you are using bioconvert, use --force ")
                raise ValueError(msg)
            else:
                _log.warning("output file will be overwritten")

        # Only fastq files can be compressed with dsrc
        if outfile.endswith(".dsrc"):
            # only valid for FastQ files extension
            # dsrc accepts only .fastq file extension
            if outfile.endswith(".fastq.dsrc") is False:
                msg = "When compressing with .dsrc extension, " +\
                    "only files ending with .fastq extension are " +\
                    "accepted. This is due to the way dsrc executable +"\
                    "is implemented."
                _log.critical(msg)
                raise IOError

        # Case1: fastq.gz to fasta.bz2
        # Here, we want to decompress, convert, compress.
        # so we need the extension without .gz or .bz2
        # We should have inext set to fastq and outext
        # set to fasta.bz2
        self.inext = getext(infile, remove_compression=True)
        self.outext = getext(outfile, remove_compression=True)

        # Case 2, fastq.gz to fastq.bz2
        # data is not changed, just the type of compression, so we want
        # to keep the original extensions, here inext and outext  will contain
        # .gz and .bz2
        if self.inext == self.outext:
            _log.info("decompression/compression mode")
            self.inext = getext(infile)
            self.outext = getext(outfile)

        self.mapper = Registry()

        # From the input parameters 1 and 2, we get the module name
        try:
            if in_fmt is None:
                in_fmt = get_format_from_extension(self.inext)
            if out_fmt is None:
                out_fmt = get_format_from_extension(self.outext)
            self.in_fmt = in_fmt.upper()
            self.out_fmt = out_fmt.upper()
            _log.info("Input: %s", self.in_fmt)
            _log.info("Output: %s", self.out_fmt)
            class_converter = self.mapper[(self.in_fmt, self.out_fmt)]
            self.name = class_converter.__name__
        except KeyError:
            # This module name was not found
            # Try to find path of converters
            conv_path = self.mapper.conversion_path(self.in_fmt, self.out_fmt)
            _log.debug("path: {}".format(conv_path))
            if conv_path:
                _log.info("Direct conversion not implemented. "
                          "Chaining converters.")
                # implemented in bioconvert/core/base.py
                # using temporary files
                class_converter = make_chain([
                    (pair, self.mapper[pair]) for pair in conv_path])
            else:
                msg = "Requested input format ('%s') to output format ('%s') is not available in bioconvert" %(
                    self.in_fmt,
                    self.out_fmt,
                )
                _log.critical(msg)
                _log.critical("Use --formats to know the available formats and --help for examples")
                raise Exception(msg)

        self.converter = class_converter(infile, outfile)
        _log.info("Using {} class".format(self.converter.name))

    def __call__(self, *args, **kwargs):
        self.converter(*args, **kwargs)

    def boxplot_benchmark(self, *args, **kwargs):
        self.converter.boxplot_benchmark(*args, **kwargs)
