# -*- coding: utf-8 -*-
#
#  This file is part of Bioconvert software
#
#  Copyright (c) 2016 - Bioconvert Development Team
#
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/biokit/bioconvert
#  documentation: http://bioconvert.readthedocs.io
#
##############################################################################
import os
import time
import abc
import select
import sys
import inspect
import shutil
import subprocess
import itertools

from io import StringIO
from subprocess import Popen, PIPE

from easydev.multicore import cpu_count

import colorlog
_log = colorlog.getLogger(__name__)


from bioconvert.core.benchmark import Benchmark
from bioconvert.core.utils import generate_outfile_name
from bioconvert import logger


class ConvMeta(abc.ABCMeta):
    """
    This metaclass checks that the converter classes have

       * an attribute input_ext
       * an attribute output_ext

    This is an abstract class used by :class:`ConvBase` class.
    The standard way to build a new converter is to inherits from :class:`ConvBase`
    or a subclasses of it, for instance::

        class Fasta_2_Fasta(ConvBase):

            input_ext = ['.fa', '.fst', '.fasta']
            output_ext = '.fa'

            def __call__(self, *args, **kwargs):
                # do conversion here
                pass

    The declaration of input_ext and output_ext is quite permissive. You can add
    prefix the extension with a dot or not; if the input consists of a single
    extension, it can be a single string, or a list/set/tuple of strings.

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

        # do not check extension since modules does not require to specify extension anymore

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
                try:
                    # get the isusable method when provided
                    isusable_method = getattr(cls, "_isusable_method_{}".format(conv_meth))
                    try:
                        # evaluate the usability of the method
                        isusable = isusable_method()
                    except Exception as e:
                        _log.debug("converter '{}' failed to evaluat '{}': {}".format(cls.__name__, conv_meth, e))
                        # the isusable method failed, the associated method should not be used
                        isusable = False
                except Exception as e:
                    # no isusable method have been found, method is supposed to work
                    isusable = True
                if isusable:
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
    """
    This is the base class for all converters.
    To build a new converter create a new class which inherits of :class:`ConvBase`
    and implement __call__ method (which is abstract). The class attributes
    input_ext and output_ext must be also override in the subclass.
    for instance: ::

        class Fasta_2_Fasta(ConvBase):

            input_ext = ['.fa', '.fst', '.fasta']
            output_ext = '.fa'

        __call__(self, *args, **kwargs):
            do conversion
    """
    # specify the extensions of the input file, can be a sequence (must be
    # overridden in subclasses)
    input_ext = None

    # specify the extensions of the output file, can be a sequence (must be
    # overridden in subclasses)
    output_ext = None
    _default_method = None
    _is_compressor = False

    def __init__(self, infile, outfile):
        """.. rubric:: constructor

        :param str infile: The path of the input file.
        :param str outfile: The path of The output file
        """
        # do not check the existence of the input file because it could be just a prefix
        # if os.path.exists(infile) is False:
        #     msg = "Incorrect input file: %s" % infile
        #     _log.error(msg)
        #     raise ValueError(msg)

        if not outfile:
            outfile = generate_outfile_name(infile, self.output_ext[0])



        self.infile = infile
        self.outfile = outfile
        self.threads = cpu_count()
        # self.max_threads = cpu_count() #TODO KO/OK ?
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

        _log.info("Executing {} method".format(method_name))
        # reference to the method requested
        method_reference = getattr(self, "_method_{}".format(method_name))

        # call the method itself
        method_reference(*args, **kwargs)


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
        t1 = time.time()
        _log.info("{}> ".format(self.name))
        _log.info("CMD: {}".format(cmd))

        shell(cmd)

        t2 = time.time()
        self.last_duration = t2 - t1
        _log.info("Took {} seconds ".format(t2-t1))
        self._last_time = t2 - t1

    def execute(self, cmd, ignore_errors=False, verbose=False, shell=False):

        if shell is True or self._execute_mode == "shell":
            self.shell(cmd)
            return

        t1 = time.time()
        _log.info("{}> ".format(self.name))
        _log.info("CMD: {}".format(cmd))
        self._execute(cmd, ignore_errors, verbose)
        t2 = time.time()
        self.last_duration = t2 - t1
        _log.info("Took {} seconds ".format(t2-t1))
        self._last_time = t2 - t1

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
        else:
            return self._default_method
    default = property(_get_default_method)

    def install_tool(self, executable):
        """Install the given tool, using the script:
        bioconvert/install_script/install_executable.sh
        if the executable is not already present

        :param executable to install
        :return: nothing
        """
        import bioconvert
        from bioconvert import bioconvert_data

        if shutil.which(executable) is None:
            logger.info("Installing tool : "+executable)
            bioconvert_path = bioconvert.__path__[0]
            script = bioconvert_data('install_'+executable+'.sh', where="../misc")
            subprocess.call(['sh',script])

    @classmethod
    def add_argument_to_parser(cls, sub_parser):
        for arg in itertools.chain(cls.get_common_arguments(), cls.get_additional_arguments()):
            arg.add_to_sub_parser(sub_parser)

    @classmethod
    def get_additional_arguments(cls):
        return []

    @classmethod
    def get_common_arguments(cls):
        yield ConvArg(
            names="input_file",
            default=None,
            help="The path to the file to convert.",
        )
        yield ConvArg(
            names="output_file",
            nargs="?",
            default=None,
            help="The path where the result will be stored.",
        )
        # yield ConvArg(
        #     names=["-i", "--input-format", ],
        #     nargs="?",
        #     default=None,
        #     help="Provide the input format. Check the --formats to see valid input name",
        # )
        # yield ConvArg(
        #     names=["-o", "--output-format", ],
        #     nargs="?",
        #     default=None,
        #     help="Provide the output format. Check the --formats to see valid input name",
        # )
        try:
            #Some converter does not have any method and work in __call__, so preventing to crash by searching for them
            yield ConvArg(
                names=["-c", "--method", ],
                nargs="?",
                default=None,
                help="The method to use to do the conversion. Default method is '%s'" % cls._get_default_method(cls),
                choices=cls.available_methods,
            )
        except Exception as e:
            _log.warning("converter '{}' does not seems to have methods: {}".format(cls.__name__, e))
            pass
        yield ConvArg(
            names=["-f", "--force", ],
            action="store_true",
            help="if outfile exists, it is overwritten with this option",
        )
