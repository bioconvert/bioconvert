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
"""Main factory of Bioconvert"""
import copy
import time
import abc
import select
import sys
import inspect
import shutil
import subprocess
import itertools

from subprocess import Popen, PIPE
from io import StringIO
from collections import deque

from easydev import TempFile
from easydev.multicore import cpu_count

import colorlog

import bioconvert

from bioconvert.core.benchmark import Benchmark
from bioconvert.core import extensions

from bioconvert.core.utils import generate_outfile_name
from bioconvert import logger
from . import BioconvertError


_log = colorlog.getLogger(__name__)


class ConvMeta(abc.ABCMeta):
    """This metaclass checks that the converter classes have

       * an attribute input_ext
       * an attribute output_ext

    This is a meta class used by :class:`ConvBase` class. For developers
    only.
    """

    @classmethod
    def split_converter_to_format(cls, converter_name: str):
        converter_name = converter_name.replace("_to_", "2")
        if '2' not in converter_name:
            raise TypeError("converter's name '{}' name must follow convention input2output".format(converter_name))
        # for BZ2 2 GZ
        if "22" in converter_name:
            input_fmt, output_fmt = converter_name.upper().split('22', 1)
            input_fmt += "2"
            input_fmt = tuple([input_fmt])
            output_fmt = tuple([output_fmt])
        else:
            input_fmt, output_fmt = converter_name.upper().split('2', 1)
            input_fmt = input_fmt.upper().split("_")
            input_fmt = tuple(input_fmt)
            output_fmt = output_fmt.upper().split("_")
            output_fmt = tuple(output_fmt)

        return input_fmt, output_fmt

    @classmethod
    def lower_tuple(cls, format_tuple):
        format_tuple = [format.lower() for format in format_tuple]
        return format_tuple


    def __init__(cls, name, bases, classdict):

        # do not check extension since modules does not require to specify
        # extension anymore

        def is_conversion_method(item):
            """Return True if method name starts with _method_

            This method is used to keep methods that starts with _method_.
            It uses inspect.getmembers func to list
            all conversion methods implemented in a convertor class.

            :param item: the object to inspect
            :return: True if method's name starts with '__method_', False otherwise.
            :rtype: boolean
            """

            return inspect.isfunction(item) and \
                 item.__name__.startswith('_method_') and \
                 item.__name__ != "_method_dummy"

        if name != 'ConvBase':
            input_fmt, output_fmt = cls.split_converter_to_format(name)
            setattr(cls, 'input_fmt', input_fmt)
            setattr(cls, 'output_fmt', output_fmt)

            if not cls.input_ext:
                # We add all the extensions for each converter into a list.
                input_ext = []
                cls.input_ext = cls.lower_tuple(cls.input_fmt)
                for format in cls.input_ext:
                    input_ext.append(tuple(extensions.extensions[format]))
                # then we turn the list into tuple as output_ext attribute
                setattr(cls, 'input_ext', tuple(input_ext))
            # if the developer did not specify an output_ext attribute
            if not cls.output_ext:
                # We add all the extensions for each converter into a list.
                output_ext = []
                cls.output_ext = cls.lower_tuple(cls.output_fmt)
                for format in cls.output_ext:
                    output_ext.append(tuple(extensions.extensions[format]))
                # then we turn the list into tuple as output_ext attribute
                setattr(cls, 'output_ext', tuple(output_ext))
                # if the key is not in the dictionary return an error message
            available_conv_meth = []
            for name in inspect.getmembers(cls, is_conversion_method):
                # do not use strip() but split()
                conv_meth = name[0].split("_method_")[1]
                is_disabled = getattr(name[1], "is_disabled", None)
                if is_disabled is None:
                    _log.debug("converter '{}': method {} is not decorated, we expect it to work all time".format(
                        cls.__name__,
                        conv_meth,
                    ))
                    is_disabled = False
                if not is_disabled:
                    available_conv_meth.append(conv_meth)
                else:
                    _log.warning("converter '{}': method {} is not available".format(cls.__name__, conv_meth, ))
            setattr(cls, 'available_methods', available_conv_meth)
            _log.debug("class = {}  available_methods = {}".format(cls.__name__, available_conv_meth))


class ConvArg(object):
    """This class can be used to add specific extra arguments to any converter

    For instance, imagine a conversion named **A2B** that requires the
    user to provide a reference. Then, you may want to provide the
    `--reference` extra argument. This is possible by adding a class
    method named get_additional_arguments that will yield instance of
    this class for each extra argument.

    ::

        @classmethod
        def get_additional_arguments(cls):
            yield ConvArg(
                names="--reference",
                default=None,
                help="the referenc"
            )

    Then, when calling bioconvert as follows,::

        bioconvert A2B --help

    the new argument will be shown in the list of arguments.


    """


    black_listed_argument_for_argparse = [
        "output_argument",
    ]

    def __init__(self, names, help, **kwargs):
        if isinstance(names, list):
            self.args_for_sub_parser = names
        else:
            self.args_for_sub_parser = [names, ]
        self.kwargs_for_sub_parser = {
            'help': help
        }
        self.kwargs_for_sub_parser.update(kwargs)

    def add_to_sub_parser(self, sub_parser):
        kwargs = copy.deepcopy(self.kwargs_for_sub_parser)
        for a in self.black_listed_argument_for_argparse:
            kwargs.pop(a, None)
        sub_parser.add_argument(*self.args_for_sub_parser, **kwargs)

    @classmethod
    def file(cls, path):
        return path


class ConvBase(metaclass=ConvMeta):
    """Base class for all converters.

    To build a new converter, create a new class which inherits from
    :class:`ConvBase` and implement method that performs the conversion.
    The name of the converter method must start with ``_method_``.

    For instance: ::

        class FASTQ2FASTA(ConvBase):

            def _method_python(self, *args, **kwargs):
                # include your code here. You can use the infile and outfile
                # attributes.
                self.infile
                self.outfile

    """
    # specify the extensions of the input file, can be a sequence (must be
    # overridden in subclasses)
    input_ext = None

    # specify the extensions of the output file, can be a sequence (must be
    # overridden in subclasses)
    output_ext = None

    # list available methods
    available_methods = []

    # default method should be provided
    _default_method = None
    _is_compressor = False
    # Can be overriden and if True, new argument --thread is added automatically
    _threading = False
    _extra_arguments = ""

    # threads to be used by default if argument is required in a method
    # this will be overriden if _threading set to True and therefore --threads
    # set by the user. It is feed back into Bioconvert class
    threads = cpu_count()

    def __init__(self, infile, outfile):
        """.. rubric:: constructor

        :param str infile: the path of the input file.
        :param str outfile: the path of The output file
        """

        self.infile = infile
        self.outfile = outfile

        # execute mode can be shell or subprocess.
        self._execute_mode = "shell"

        # The logger to be set to INFO, DEBUG, WARNING, ERROR, CRITICAL
        self.logger = logger

    def __call__(self, *args, method_name=None, **kwargs):
        """

        :param str method_name: the method to be found in :attr:`available_methods`
        :param str method: the method to be found in :attr:`available_methods`
        :param *args: positional arguments
        :param *kwargs: keyword arguments

        """
        # If method provided, use it
        if "method" in kwargs:
            method_name = kwargs["method"]
            del kwargs["method"]

        # If not, but there is one argument, presumably this is
        # the method
        method_name = method_name or self.default

        # If not, we need to check the name
        # "dummy" is a method used to evaluate the cost of the
        # execute() method for the benchmark
        if method_name not in self.available_methods + ['dummy']:
            msg = "Methods available are {}".format(self.available_methods)
            _log.error(msg)
            raise ValueError(msg)

        _log.info("{}> Executing {} method ".format(self.name, method_name))
        # reference to the method requested
        method_reference = getattr(self, "_method_{}".format(method_name))

        # call the method itself

        t1 = time.time()
        method_reference(*args, **kwargs)
        t2 = time.time()
        _log.info("Took {} seconds ".format(t2 - t1))

    @property
    def name(self):
        """
        The name of the class
        """
        return type(self).__name__

    def _method_dummy(self, *args, **kwargs):
        # The execute commands has a large initialisation cost (about a second)
        # This commands does not and can be used to evaluate that cost
        self.execute("")

    def shell(self, cmd):
        from bioconvert.core.shell import shell
        _log.info("CMD: {}".format(cmd))
        shell(cmd)

    def execute(self, cmd, ignore_errors=False, verbose=False, shell=False):

        if ">" in cmd:
            lhs, rhs = cmd.split(">", 1)
            cmd = lhs + self._extra_arguments + ">" + rhs
        else:
            cmd = cmd + self._extra_arguments

        if shell is True or self._execute_mode == "shell":
            self.shell(cmd)
            return
        self._execute(cmd, ignore_errors, verbose)

    def _execute(self, cmd, ignore_errors=False, verbose=False):
        """
        Execute a command in a sub-shell

        :param str cmd: the command to execute
        :param ignore_errors: If True the result is returned whatever the
                              return value of the sub-shell.
                              Otherwise a Runtime error is raised when the sub-shell
                              return a non zero value
        :param verbose: If true displays errors on standard error
        :return: the result of the command
        :rtype: a :class:`StringIO` instance
        """
        try:
            process_ = Popen(cmd,
                             shell=True,
                             stdout=PIPE,
                             stderr=PIPE,
                             stdin=None)
        except Exception as err:
            msg = "Failed to execute Command: '{}'. error: '{}'".format(cmd, err)
            raise RuntimeError(msg)

        inputs = [process_.stdout, process_.stderr]
        output = StringIO()
        errors = StringIO()
        while process_.poll() is None:
            # select has 3 parameters, 3 lists, the sockets, the fileobject to watch
            # in reading, writing, the errors
            # in addition a timeout option (the call is blocking while a fileObject
            # is not ready to be processed)
            # by return we get 3 lists with the fileObject to be processed
            # in reading, writing, errors.
            readable, writable, exceptional = select.select(inputs, [], [], 1)

            while readable and inputs:
                for flow in readable:
                    data = flow.read()
                    if not data:
                        # the flow ready in reading which has no data
                        # is a closed flow
                        # thus we must stop to watch it
                        inputs.remove(flow)
                    if flow is process_.stdout:
                        output.write(data.decode("utf-8"))
                    elif flow is process_.stderr:
                        errors.write(data.decode("utf-8"))
                        print(process_.stderr)
                readable, writable, exceptional = select.select(inputs, [], [], 1)

        errors = errors.getvalue().strip()
        if verbose:
            if errors:
                print(errors, file=sys.stderr)

        if process_.returncode != 0:
            if not ignore_errors:
                raise RuntimeError(errors)
        else:
            return output

    def boxplot_benchmark(self, N=5, rerun=True, include_dummy=False,
                          to_exclude=[], to_include=[], rot_xticks=90,
                          boxplot_args={}):
        """Simple wrapper to call :class:`Benchmark` and plot the results

        see :class:`~bioconvert.core.benchmark.Benchmark` for details.

        """
        if to_include == "all":
            to_include = []

        self._benchmark = Benchmark(self, N=N, to_exclude=to_exclude,
                                    to_include=to_include)
        self._benchmark.include_dummy = include_dummy
        data = self._benchmark.plot(rerun=rerun, rot_xticks=rot_xticks,
                                    boxplot_args=boxplot_args)
        return data

    def _get_default_method(self):
        if self._default_method is None:
            return self.available_methods[0]
        elif self._default_method not in self.available_methods:
            return self.available_methods[0]
        else:
            return self._default_method
    default = property(_get_default_method)

    def install_tool(self, executable):
        """Install the given tool, using the script:
        bioconvert/install_script/install_executable.sh
        if the executable is not already present

        :param executable: executable to install
        :return: nothing

        """
        # imported but not unused (when we don't have bioconvert_path)
        # import bioconvert
        from bioconvert import bioconvert_data

        if shutil.which(executable) is None:
            logger.info("Installing tool : " + executable)
            # Assigned but never used, says flake8
            # bioconvert_path = bioconvert.__path__[0]
            script = bioconvert_data(
                'install_' + executable + '.sh', where="../misc")
            subprocess.call(['sh', script])

    @classmethod
    def add_argument_to_parser(cls, sub_parser):
        sub_parser.description = cls.get_description()
        for arg in itertools.chain(cls.get_IO_arguments(),
                                   cls.get_common_arguments_for_converter(),
                                   cls.get_additional_arguments()):
            arg.add_to_sub_parser(sub_parser)

    @classmethod
    def get_description(cls):
        msg = "Convert file from '{}' to '{}' format. "
        msg += "See bioconvert.readthedocs.io for details" 
        msg = msg.format(*ConvMeta.split_converter_to_format(cls.__name__))
        return msg

    @classmethod
    def get_additional_arguments(cls):
        return []


    # common arguments for the sub command case
    # when using bioconvert <conversion>
    @staticmethod
    def get_IO_arguments():
        yield ConvArg(
            names="input_file",
            nargs="?",
            default=None,
            type=ConvArg.file,
            help="The path to the file to convert.",
        )
        yield ConvArg(
            names="output_file",
            nargs="?",
            default=None,
            type=ConvArg.file,
            output_argument=True,
            help="The path where the result will be stored.",
        )

    @staticmethod
    def get_common_arguments():
        yield ConvArg(
            names=["-f", "--force", ],
            action="store_true",
            help="if outfile exists, it is overwritten with this option",
        )
        yield ConvArg(
            names=["-v", "--verbosity", ],
            default=bioconvert.logger.level,
            help="Set the outpout verbosity.",
            choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        )
        yield ConvArg(
            names=["--raise-exception", ],
            action="store_true",
            help="Let exception ending the execution be raised and displayed",
        )
        yield ConvArg(
            names=["-X", "--batch", ],
            default=False,
            action="store_true",
            help="Allow conversion of a set of files using wildcards. You "
                 "must use quotes to escape the wildcards. For instance: "
                 "--batch 'test*fastq' ")
        yield ConvArg(
            names=["-b", "--benchmark", ],
            default=False,
            action="store_true",
            help="Running all available methods",
        )
        yield ConvArg(
            names=["-N", "--benchmark-N", ],
            default=5,
            type=int,
            help="Number of trials for each methods",
        )
        yield ConvArg(
            names=["-B", "--benchmark-methods", ],
            default="all",
            nargs="+",
            type=str,
            help="Methods to include",
        )
        yield ConvArg(
            names=["-a", "--allow-indirect-conversion", ],
            default=False,
            action="store_true",
            help="Allow to chain converter when direct conversion is absent",
        )
        yield ConvArg(
            names=["-e", "--extra-arguments", ],
            default="",
            help="Any arguments accepted by the method's tool",
        )


    @classmethod
    def get_common_arguments_for_converter(cls):
        for a in ConvBase.get_common_arguments():
            yield a
        try:
            # Some converters do not have any method and work
            # in __call__, so preventing to crash by searching for them
            yield ConvArg(
                names=["-m", "--method", ],
                nargs="?",
                default=cls._get_default_method(cls),
                help="The method to use to do the conversion.",
                choices=cls.available_methods,
            )
        except Exception as e:
            _log.warning("converter '{}' does not seems to have methods: {}".format(cls.__name__, e))
            pass
        yield ConvArg(
            names=["-s", "--show-methods", ],
            default=False,
            action="store_true",
            help="A converter may have several methods",
        )

        if cls._threading:
            yield ConvArg(
               names=["-t", "--threads"],
               #nargs=1,
               type=int,
               default=cls.threads,
               help="threads to be used",
            )


# Implementing a class creator
# The created class will have the correct name, will inherit from ConvBase
# It will have a conversion method chaining conversions through tempfiles
def make_chain(converter_map):
    """
    Create a class performing step-by-step conversions following a path.
    *converter_map* is a list of pairs ((in_fmt, out_fmt), converter).
    It describes the conversion path.
    """
    in_fmt = converter_map[0][0][0]
    out_fmt = converter_map[-1][0][1]
    chain_name = "{}2{}".format("_".join(in_fmt), "_".join(out_fmt))
    chain_attributes = {}

    def chain_init(self, infile, outfile):
        super().__init__(infile, outfile)
        self._default_method = "chain"

    def _method_chain(self, *args, **kwargs):
        """This method successively uses the default conversion method of each
        converter in the conversion path."""

        def conv_step(converter, infile, outfile):
            """Performs one conversion step."""
            converter(infile, outfile)(*args, **kwargs)

        # Contains the last temporary output file, if any
        pipe_files = deque()
        for (step_num, ((_, out_fmt), converter)) \
                in enumerate(self.converter_map, start=1):
            if step_num == 1:
                # May not be necessary:
                step_infile = None
                step_input = self.infile
                del_infile = False
            else:
                step_infile = pipe_files.popleft()
                step_input = step_infile.name
                del_infile = True

            if step_num == self.nb_steps:
                # May not be necessary:
                step_outfile = None
                step_output = self.outfile
            else:

                #FIXME: for mutiple IO converters
                if len(out_fmt) == 1:
                    step_outfile = TempFile(suffix=out_fmt[0].lower())
                    step_output = step_outfile.name
                    pipe_files.append(step_outfile)

            conv_step(converter, step_input, step_output)
            if del_infile:
                step_infile.delete()

    chain_attributes["converter_map"] = converter_map
    chain_attributes["nb_steps"] = len(converter_map)
    chain_attributes["__init__"] = chain_init
    chain_attributes["_method_chain"] = _method_chain
    chain = type(chain_name, (ConvBase,), chain_attributes)
    # https://stackoverflow.com/a/43779009/1878788
    # Allows calling super in chain.__init__
    __class__ = chain
    return chain
