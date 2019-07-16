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
"""Provides a general tool to perform pre/post compression"""
from distutils.spawn import find_executable
from functools import wraps
from os.path import splitext

import colorlog
import pkg_resources
from easydev import TempFile

_log = colorlog.getLogger(__name__)


def in_gz(func):
    """Marks a function as accepting gzipped input."""
    func.in_gz = True
    return func


def make_in_gz_tester(converter):
    """Generates a function testing whether a conversion method of *converter*
    has the *in_gz* tag."""

    def is_in_gz(method):
        """Accesses the function corresponding to *method* and tells whether it
        has the *in_gz* tag."""
        return hasattr(getattr(
            converter, "_method_{}".format(method)), "in_gz")

    return is_in_gz


def compressor(func):
    """Decompress/compress input file without pipes

    Does not use pipe: we decompress and compress back the input file.
    The advantage is that it should work for any files (even very large).

    This decorator should be used by method that uses pure python code
    """

    # https://stackoverflow.com/a/309000/1878788
    @wraps(func)
    def wrapped(inst, *args, **kwargs):
        infile_name = inst.infile
        if type(inst.outfile) is not list:

            output_compressed = None
            if inst.outfile.endswith(".gz"):
                (inst.outfile, output_compressed) = splitext(inst.outfile)
            elif inst.outfile.endswith(".bz2"):
                (inst.outfile, output_compressed) = splitext(inst.outfile)
            elif inst.outfile.endswith(".dsrc"):  # !!! only for fastq files
                (inst.outfile, output_compressed) = splitext(inst.outfile)
            # Now inst has the uncompressed output file name

        if infile_name.endswith(".gz"):
            # decompress input
            # TODO: https://stackoverflow.com/a/29371584/1878788
            _log.info("Generating uncompressed version of {} ".format(infile_name))
            (ungz_name, _) = splitext(infile_name)
            (_, base_suffix) = splitext(ungz_name)
            with TempFile(suffix=base_suffix) as ungz_infile:
                inst.infile = ungz_infile.name
                inst.shell("unpigz -c -p {} {} > {}".format(
                    inst.threads, infile_name, inst.infile))
                # computation
                results = func(inst, *args, **kwargs)
            inst.infile = infile_name
        else:
            results = func(inst, *args, **kwargs)

        # Compress output and restore inst output file name
        if output_compressed == ".gz":
            # TODO: this uses -f ; should be a
            _log.info("Compressing output into .gz")
            inst.shell("pigz -f -p {} {}".format(inst.threads, inst.outfile))
            inst.outfile = inst.outfile + ".gz"
        elif output_compressed == ".bz2":
            _log.info("Compressing output into .bz2")
            inst.shell("pbzip2 -f -p{} {}".format(inst.threads, inst.outfile))
            inst.outfile = inst.outfile + ".bz2"
        elif output_compressed == ".dsrc":  # !!! only for FastQ files
            _log.info("Compressing output into .dsrc")
            inst.shell("dsrc c -t{} {} {}.dsrc".format(
                inst.threads, inst.outfile, inst.outfile))
            inst.outfile = inst.outfile + ".dsrc"
        return results

    return in_gz(wrapped)


def out_compressor(func):
    """Compress output file without pipes

    This decorator should be used by method that uses pure python code
    """

    # https://stackoverflow.com/a/309000/1878788
    @wraps(func)
    def wrapped(inst, *args, **kwargs):
        output_compressed = None
        if inst.outfile.endswith(".gz"):
            (inst.outfile, output_compressed) = splitext(inst.outfile)
        elif inst.outfile.endswith(".bz2"):
            (inst.outfile, output_compressed) = splitext(inst.outfile)
        elif inst.outfile.endswith(".dsrc"):  # !!! only for fastq files
            (inst.outfile, output_compressed) = splitext(inst.outfile)
        # Now inst has the uncompressed output file name

        # computation
        results = func(inst, *args, **kwargs)

        # Compress output and restore inst output file name
        if output_compressed == ".gz":
            # TODO: this uses -f ; should be a
            _log.info("Compressing output into .gz")
            inst.shell("pigz -f -p {} {}".format(inst.threads, inst.outfile))
            inst.outfile = inst.outfile + ".gz"
        elif output_compressed == ".bz2":
            _log.info("Compressing output into .bz2")
            inst.shell("pbzip2 -f -p{} {}".format(inst.threads, inst.outfile))
            inst.outfile = inst.outfile + ".bz2"
        elif output_compressed == ".dsrc":  # !!! only for FastQ files
            _log.info("Compressing output into .dsrc")
            inst.shell("dsrc c -t{} {} {}.dsrc".format(
                inst.threads, inst.outfile, inst.outfile))
            inst.outfile = inst.outfile + ".dsrc"
        return results

    return wrapped


def requires_nothing(func):
    """Marks a function as not needing dependencies."""
    func.is_disabled = False
    return func


def requires(
        external_binary=None,
        python_library=None,
        external_binaries=None,
        python_libraries=None,
):
    """

    :param external_binary: a system binary required for the method
    :param python_library:  a python library required for the method
    :param external_binaries: an array of system binaries required for the method
    :param python_libraries: an array of python libraries required for the method
    :return:
    """
    external_binaries = external_binaries or []
    python_libraries = python_libraries or []
    if external_binary:
        external_binaries.append(external_binary)
    if python_library:
        python_libraries.append(python_library)

    __missing_binaries = getattr(requires, "__missing_binaries", {})
    requires.__missing_binaries = __missing_binaries
    __missing_libraries = getattr(requires, "__missing_libraries", {})
    requires.__missing_libraries = __missing_libraries
    __pip_libraries = getattr(requires, "__pip_libraries", None)
    if __pip_libraries is None:
        __pip_libraries = [p.project_name for p in pkg_resources.working_set]
        requires.__pip_libraries = __pip_libraries

    def real_decorator(function):
        @wraps(function)
        def wrapped(inst, *args, **kwargs):
            return function(inst, *args, **kwargs)

        try:
            for bin in external_binaries:
                try:
                    if __missing_binaries[bin]:
                        raise Exception("{} has already be seen as missing".format(bin))
                except KeyError:
                    __missing_binaries[bin] = True
                    # shell("which %s" % bin)
                    if find_executable(bin) is None:
                        raise Exception("{} was not found in path".format(bin))
                    __missing_binaries[bin] = False
            for lib in python_libraries:
                try:
                    if __missing_libraries[lib]:
                        raise Exception("{} has already be seen as missing".format(lib))
                except KeyError:
                    missing = lib not in __pip_libraries
                    __missing_libraries[lib] = missing
                    if missing:
                        raise Exception("{} was not found by pip".format(lib))
            wrapped.is_disabled = False
        except Exception as e:
            _log.debug(e)
            wrapped.is_disabled = True
        return wrapped

    return real_decorator


def get_known_dependencies_with_availability(as_dict=False):
    if as_dict:
        external_binaries = {}
        python_libraries = {}
        for binary, missing in getattr(requires, "__missing_binaries", {}).items():
            external_binaries[binary] = dict(
                available=not missing,
            )
        for library, missing in getattr(requires, "__missing_libraries", {}).items():
            python_libraries[library] = dict(
                available=not missing,
            )
        return dict(
            external_binaries=external_binaries,
            python_libraries=python_libraries,
        )
    ret = []
    for binary, status in sorted(getattr(requires, "__missing_binaries", {}).items()):
        ret.append((binary, not status, "binary",))
    for library, status in sorted(getattr(requires, "__missing_libraries", {}).items()):
        ret.append((library, not status, "library",))
    return ret
