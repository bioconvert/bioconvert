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
"""Provides a general tool to perform pre/post compression"""
from functools import wraps
from bioconvert import logger


def compressor(func):
    """Decompress/compress input file without pipes 

    Does not use pipe: we decompress and compress back the input file.
    The advantage is that it should work for any files (even very large).

    This decorator should be used by method that uses pure python code
    """
    @wraps(func)
    def wrapped(inst, *args, **kwargs):
        infile_name = inst.infile

        output_compressed = None
        if inst.outfile.endswith(".gz"):
            output_compressed = ".gz"
            inst.outfile = inst.outfile.split(".gz")[0]
        elif inst.outfile.endswith(".bz2"):
            output_compressed = ".bz2"
            inst.outfile = inst.outfile.split(".bz2")[0]
        elif inst.outfile.endswith(".dsrc"): # !!! only for fastq files
            output_compressed = ".dsrc"
            inst.outfile = inst.outfile.split(".dsrc")[0]

        if inst.infile.endswith(".gz"):
            # decompress input
            logger.info("Decompressing %s " % inst.infile)
            newinfile = inst.infile.split(".gz")[0]
            inst.shell("unpigz -p {} {}".format(inst.threads, inst.infile))
            inst.infile = newinfile
            # computation
            results = func(inst, *args, **kwargs)
            logger.info("Compressing back %s" % inst.infile)
            # compress back the input
            inst.shell("pigz -p {} {}".format(inst.threads, inst.infile))
            inst.infile = infile_name
        else:
            results = func(inst, *args, **kwargs)

        if output_compressed == ".gz":
            # TODO: this uses -f ; should be a
            logger.info("Compressing output into .gz")
            inst.shell("pigz -f -p {} {}".format(inst.threads, inst.outfile))
            inst.outfile = inst.outfile + ".gz"
        elif output_compressed == ".bz2":
            logger.info("Compressing output into .bz2")
            inst.shell("pbzip2 -f -p{} {}".format(inst.threads, inst.outfile))
            inst.outfile = inst.outfile + ".bz2"
        elif output_compressed == ".dsrc":   # !!! only for FastQ files
            logger.info("Compressing output into .dsrc")
            inst.shell("dsrc c  -t{} {} {}.dsrc".format(inst.threads,
                inst.outfile, inst.outfile))
            inst.outfile = inst.outfile + ".dsrc"
        return results
    return wrapped

