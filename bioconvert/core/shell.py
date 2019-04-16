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
"""Simplified version of shell.py module from snakemake package

Temporary module that would probably need some checking / simplification

This is a temporary replacement to the execute() in convbase
that also uses subprocess.Popen but is slower
"""
import _io
import sys
import os
import subprocess as sp

import colorlog
_log = colorlog.getLogger(__name__)


STDOUT = sys.stdout


class shell:
    _process_args = {}
    _process_prefix = ""
    _process_suffix = ""

    @classmethod
    def executable(cls, cmd):
        if os.path.split(cmd)[-1] == "bash":
            cls._process_prefix = "set -euo pipefail; "
        cls._process_args["executable"] = cmd

    @classmethod
    def prefix(cls, prefix):
        cls._process_prefix = prefix 

    @classmethod
    def suffix(cls, suffix):
        cls._process_suffix = suffix

    def __new__(cls, cmd, *args,
                is_async=False,
                iterable=False,
                read=False, **kwargs):
        if "stepout" in kwargs:
            raise KeyError("Argument stepout is not allowed in shell command.")

        stdout = sp.PIPE if iterable or is_async or read else STDOUT

        close_fds = sys.platform != 'win32'

        _log.debug(cmd)
        proc = sp.Popen("{} {} {}".format(
                            cls._process_prefix,
                            cmd.rstrip(),
                            cls._process_suffix),
                        bufsize=-1,
                        shell=True,
                        stdout=stdout,
                        close_fds=close_fds, **cls._process_args)

        ret = None
        if iterable:
            return cls.iter_stdout(proc, cmd)
        if read:
            ret = proc.stdout.read()
        elif is_async:
            return proc
        retcode = proc.wait()
        if retcode:
            raise sp.CalledProcessError(retcode, cmd)
        return ret

    @staticmethod
    def iter_stdout(proc, cmd):
        for l in proc.stdout:
            yield l[:-1].decode()
        retcode = proc.wait()
        if retcode:
            raise sp.CalledProcessError(retcode, cmd)
