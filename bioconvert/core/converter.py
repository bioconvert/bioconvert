# -*- coding: utf-8 -*-

###########################################################################
# Bioconvert is a project to facilitate the interconversion               #
# of life science data from one format to another.                        #
#                                                                         #
# Authors: see CONTRIBUTORS.rst                                           #
# Copyright © 2018  Institut Pasteur, Paris and CNRS.                     #
# See the COPYRIGHT file for details                                      #
#                                                                         #
# bioconvert is free software: you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or       #
# (at your option) any later version.                                     #
#                                                                         #
# bioconvert is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of          #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
# GNU General Public License for more details.                            #
#                                                                         #
# You should have received a copy of the GNU General Public License       #
# along with this program (COPYING file).                                 #
# If not, see <http://www.gnu.org/licenses/>.                             #
###########################################################################
""".. rubric:: Standalone application dedicated to conversion"""
import os

import colorlog

from bioconvert.core.base import ConvMeta
from bioconvert.core.registry import Registry

_log = colorlog.getLogger(__name__)

from bioconvert.core.base import make_chain
from bioconvert.core.utils import get_extension as getext
from bioconvert.core.utils import get_format_from_extension
import sys

__all__ = ['Bioconvert']


class Bioconvert(object):
    """Universal converter used by the standalone

    ::

        from bioconvert import Bioconvert
        c = Bioconvert("test.fastq", "test.fasta", threads=4, force=True)


    """
    def __init__(self, infile, outfile, force=False,
            threads=None, extra=None):
        """.. rubric:: constructor

        :param str infile: The path of the input file.
        :param str outfile: The path of The output file
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

        if type(outfile) is str:
            outfile = [outfile]

        if type(infile) is str:
            infile = [infile]

        # some checking on the output files (existence, special case of dsrc)
        for filename in outfile:
            if os.path.exists(filename) is True:
                msg = "output file {} exists already.".format(filename)
                if force is False:
                    _log.critical("output file exists. If you are using bioconvert, use --force ")
                    raise ValueError(msg)
                else:
                    _log.warning(msg + " --force used so will be over written")

            # Only fastq files can be compressed with dsrc
            if filename.endswith(".dsrc"):
                # only valid for FastQ files extension
                # dsrc accepts only .fastq file extension
                if filename.endswith(".fastq.dsrc") is False:
                    msg = "When compressing with .dsrc extension, " +\
                        "only files ending with .fastq extension are " +\
                        "accepted. This is due to the way dsrc executable "+\
                        "is implemented."
                    _log.critical(msg)
                    raise IOError

        Lin = len(infile)
        Lout = len(outfile)

        self.inext = []
        self.outext = []
        # populate the inext
        for filename in infile:
            # example: fastq.gz to fasta.bz2
            # Here, we want to decompress, convert, compress.
            # so we need the extension without .gz or .bz2
            # We should have inext set to fastq and outext
            # set to fasta.bz2
            self.inext.append(getext(filename, remove_compression=True))

        # populate the outext 
        for filename in outfile:
            self.outext.append(getext(filename, remove_compression=True))

        # special case one to one for compression/decompression
        # Case 2, fastq.gz to fastq.bz2
        # data is not changed, just the type of compression, so we want
        # to keep the original extensions, here inext and outext  will contain
        # .gz and .bz2
        # if 1 to 1 and same extension, we overwrite self.inext and self.outext
        if Lin == Lout == 1:
            if self.inext == self.outext:
                _log.info("decompression/compression mode")
                self.inext = [getext(infile[0])]
                self.outext = [getext(outfile[0])]

        self.mapper = Registry()

        # From the input parameters 1 and 2, we get the module name
        if not list(set(list(self.mapper.get_converters_names())).intersection(sys.argv)):
            # get format from extensions
            in_fmt = [get_format_from_extension(x) for x in self.inext]
            out_fmt = [get_format_from_extension(x) for x in self.outext]
        else:
            in_fmt, out_fmt = ConvMeta.split_converter_to_format(
                list(set(list(self.mapper.get_converters_names())).intersection(sys.argv))[0])

        self.in_fmt = in_fmt
        self.out_fmt = out_fmt

        self.in_fmt = [format.lower() for format in in_fmt]
        self.in_fmt = tuple(in_fmt)

        self.out_fmt = [format.lower() for format in out_fmt]
        self.out_fmt = tuple(out_fmt)

        _log.info("Input: {}".format(self.in_fmt))
        _log.info("Output: {}".format(self.out_fmt))

        try:
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
                msg = "Requested input format ('{}') to output format ('{}') is not available in bioconvert".format(
                    self.in_fmt,
                    self.out_fmt,
                )
                _log.critical(msg)
                _log.critical("Use --formats to know the available formats and --help for examples")
                raise Exception(msg)

        # If --threads provided, we update the threads attribute


        #FIXME: hack for the compression/decompression decorators

        if Lin == 1:
            infile = infile[0]

        if Lout == 1:
            outfile = outfile[0]

        self.converter = class_converter(infile, outfile)
        if threads is not None:
            self.converter.threads = threads
        if extra:
            self.converter._extra_arguments = extra

        _log.info("Using {} class (with {} threads if needed)".format(
            self.converter.name,
            self.converter.threads))

    def __call__(self, *args, **kwargs):
        self.converter(*args, **kwargs)

    def boxplot_benchmark(self, *args, **kwargs):
        self.converter.boxplot_benchmark(*args, **kwargs)
