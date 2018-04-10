# -*- coding: utf-8 -*-
#
#  This file is part of Bioconvert software
#
#  Copyright (c) 2018 - Bioconvert Development Team
#
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/biokit/bioconvert
#  documentation: http://bioconvert.readthedocs.io
#
##############################################################################
# import os
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

_log = colorlog.getLogger(__name__)

from bioconvert.core.benchmark import Benchmark
from bioconvert.core.utils import generate_outfile_name
from bioconvert import logger


class ConvMeta(abc.ABCMeta):
    """
    This metaclass checks that the converter classes have

       * an attribute input_ext
       * an attribute output_ext

    This is an abstract class used by :class:`ConvBase` class. For developers
    only.
    """

    @classmethod
    def split_converter_to_extensions(cls, converter_name: str):
        converter_name = converter_name.replace("_to_", "2")
        if '2' not in converter_name:
            raise TypeError("converter's name '%s' name must follow convention input2output" % converter_name)
        # for BZ2 2 GZ
        if "22" in converter_name:
            input_fmt, output_fmt = converter_name.upper().split('22', 1)
            input_fmt += "2"
        else:
            input_fmt, output_fmt = converter_name.upper().split('2', 1)
        return input_fmt, output_fmt

    def __init__(cls, name, bases, classdict):

        # do not check extension since modules does not require to specify
        # extension anymore

        # def check_ext(ext, io_name):
        #     """
        #     Check if the extension is specified correctly.
        #     I must be a string or a sequence of string, otherwise raise an error
        #     it should start with a dot. Otherwise fix extension and inject it in the class

        #     :param ext: the value of the class attribute (input|output)_ext
        #     :type ext: a string or a list, tuple or set of strings
        #     :param str io_name: the type of extension, 'input' or output'
        #     :raise TypeError:  if ext is neither a string nor a sequence of strings
        #     """
        #     if isinstance(ext, str):
        #         if not ext.startswith('.'):
        #             ext = '.' + ext
        #         setattr(cls, '{}_ext'.format(io_name),  (ext, ))
        #     elif isinstance(ext, (list, tuple, set)):
        #         if not all([isinstance(one_ext, str) for one_ext in ext]):
        #             raise TypeError("each element of the class attribute '{}.{}_ext' "
        #                             "must be a string".format(cls, io_name))
        #         else:
        #             if not all([one_ext.startswith('.') for one_ext in ext]):
        #                 fixed_ext = []
        #                 for one_ext in ext:
        #                     if one_ext.startswith('.'):
        #                         fixed_ext.append(one_ext)
        #                     else:
        #                         fixed_ext.append('.' + one_ext)
        #                 setattr(cls, '{}_ext'.format(io_name), fixed_ext)
        #     else:
        #         import sys
        #         err = "the class attribute '{}.{}_ext' " \
        #               "must be specified in the class or subclasses".format(cls.__name__, io_name)
        #         _log.warning("skip class '{}': {}".format(cls.__name__, err, file=sys.stderr))
        #         raise TypeError("the class attribute '{}.{}_ext' must be specified "
        #                         "in the class or subclasses".format(cls.__name__, io_name))
        #     return True

        def is_conversion_method(item):
            """Return True is method name starts with _method_

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
            input_fmt, output_fmt = cls.split_converter_to_extensions(name.upper())
            # modules have no more input_ext and output_ext attributes
            # input_ext = getattr(cls, 'input_ext')
            # if check_ext(input_ext, 'input'):
            #     output_ext = getattr(cls, 'output_ext')
            #     check_ext(output_ext, 'output')
            setattr(cls, 'input_fmt', input_fmt)
            setattr(cls, 'output_fmt', output_fmt)
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

    def __init__(self, names, help, **kwargs):
        if isinstance(names, list):
            self.args_for_sub_parser = names
        else:
            self.args_for_sub_parser = [names, ]
        self.kwargs_for_sub_parser = dict(
            help=help,
            **kwargs,
        )

    def add_to_sub_parser(self, sub_parser):
        sub_parser.add_argument(*self.args_for_sub_parser, **self.kwargs_for_sub_parser)


class ConvBase(metaclass=ConvMeta):
    """ base class for all converters.

    To build a new converter create a new class which inherits from
    :class:`ConvBase` and implement method that performs the conversion.
    The name of the converter method must start with _method_.

    For instance: ::

        class Fastq2Fasta(ConvBase):

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
    _default_method = None
    _is_compressor = False
    threads = cpu_count()

    def __init__(self, infile, outfile):
        """.. rubric:: constructor

        :param str infile: the path of the input file.
        :param str outfile: the path of The output file
        """
        if not outfile:
            outfile = generate_outfile_name(infile, self.output_ext[0])

        self.infile = infile
        self.outfile = outfile
        self._execute_mode = "shell" #"subprocess"  # set to shell to call shell() method
        self.logger = logger

    def __call__(self, *args, method_name=None, **kwargs):
        """

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
        _log.info("Took {} seconds ".format(t2-t1))


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

        if shell is True or self._execute_mode == "shell":
            self.shell(cmd)
            return
        _log.info("CMD: {}".format(cmd))
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
            to_exclude=[], to_include=[]):
        """Simple wrapper to call :class:`Benchmark` and plot the results

        see :class:`~bioconvert.core.benchmark.Benchmark` for details.

        """
        self._benchmark = Benchmark(self, N=N, to_exclude=to_exclude,
                                    to_include=to_include)
        self._benchmark.include_dummy = include_dummy
        self._benchmark.plot(rerun=rerun)

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
        for arg in itertools.chain(cls.get_common_arguments_for_converter(), cls.get_additional_arguments()):
            arg.add_to_sub_parser(sub_parser)

    @classmethod
    def get_description(cls):
        return "Allow to convert file in '%s' to '%s' format." % ConvMeta.split_converter_to_extensions(cls.__name__)

    @classmethod
    def get_additional_arguments(cls):
        return []

    @staticmethod
    def get_common_arguments():
        yield ConvArg(
            names="input_file",
            nargs="?",
            default=None,
            help="The path to the file to convert.",
        )
        yield ConvArg(
            names="output_file",
            nargs="?",
            default=None,
            help="The path where the result will be stored.",
        )
        yield ConvArg(
            names=["-f", "--force", ],
            action="store_true",
            help="if outfile exists, it is overwritten with this option",
        )
        yield ConvArg(
            names=["-v", "--verbosity", ],
            default=bioconvert.logger.level,
            help="Set the outpout verbosity. Should be one of DEBUG, INFO, WARNING, ERROR, CRITICAL",
        )
        yield ConvArg(
            names=["--raise-exception", ],
            action="store_true",
            help="Let exception ending the execution be raised and displayed",
        )
        yield ConvArg(
            names=["-m", "--batch", ],
            default=False, action="store_true",
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
            names=["-a", "--allow-indirect-conversion", ],
            default=False,
            action="store_true",
            help="Allow to chain converter when direct conversion is absent",
        )

    @classmethod
    def get_common_arguments_for_converter(cls):
        for a in ConvBase.get_common_arguments():
            yield a
        try:
            # Some converter does not have any method and work in __call__, so preventing to crash by searching for them
            yield ConvArg(
                names=["-c", "--method", ],
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
    chain_name = "2".join([in_fmt.capitalize(), out_fmt.capitalize()])
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
                step_outfile = TempFile(suffix=out_fmt.lower())
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
