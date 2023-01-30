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
"""misc utility functions """
import os
import sys
import tempfile
import hashlib
import time

import bioconvert
from bioconvert.core.extensions import extensions

__all__ = ["get_extension", "get_format_from_extension", "generate_outfile_name", "md5", "TempFile", "Timer"]


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
        return os.path.splitext(filename)[-1].lstrip(".")


def generate_outfile_name(infile, out_extension):
    """simple utility to replace the file extension with the given one.

    :param str infile: the path to the Input file
    :param str out_extension: Desired extension
    :return: The file path with the given extension
    :rtype: str
    """
    return "{}.{}".format(os.path.splitext(infile)[0], out_extension)


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
    bioconvert.logger.critical("No format was found for extension '{}'".format(extension))
    bioconvert.logger.critical("Use --formats to know the available formats and --help" " for examples")
    sys.exit(1)


def compressor(infile, comp_ext, threads=4):
    # FIXME: could be a method in ConvBase

    if comp_ext == ".gz":
        _log.info("Compressing output into .gz")
        shell("pigz -f -p {} {}".format(threads, infile))
    elif comp_ext == ".bz2":
        _log.info("Compressing output into .bz2")
        shell("pbzip2 -f -p{} {}".format(threads, infile))
    elif comp_ext == ".dsrc":  # !!! only for FastQ files
        _log.info("Compressing output into .dsrc")
        shell("dsrc c -t{} {} {}.dsrc".format(inst.threads, inst.outfile, inst.infile))


def md5(fname, chunk=65536):
    """Return the MD5 checksums of a file

    """
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for this in iter(lambda: f.read(chunk), b""):
            hash_md5.update(this)
    return hash_md5.hexdigest()


class TempFile(object):
    """A small wrapper around tempfile.NamedTemporaryFile function

    ::

        f = TempFile(suffix="csv")
        f.name
        f.delete() # alias to delete=False and close() calls


    Copy from easydev package
    """

    def __init__(self, suffix="", dir=None):
        self.temp = tempfile.NamedTemporaryFile(suffix=suffix, delete=False, dir=dir)

    def delete(self):
        try:
            self.temp._closer.delete = True
        except:  # pragma: no cover
            self.temp.delete = True
        self.temp.close()

    def _get_name(self):
        return self.temp.name

    name = property(_get_name)

    def __exit__(self, type, value, traceback):
        try:
            self.delete()
        except AttributeError:  # pragma: no cover
            pass
        finally:
            self.delete()

    def __enter__(self):
        return self


class Timer(object):  # pragma: no cover
    """Timer working with *with* statement

    Copy from easydev package.
    """

    def __init__(self, times):
        self.times = times

    def __enter__(self):
        self.t1 = time.time()

    def __exit__(self, type, value, traceback):
        self.t2 = time.time()
        self.times.append(self.t2 - self.t1)

