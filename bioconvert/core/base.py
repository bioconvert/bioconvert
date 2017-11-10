# -*- coding: utf-8 -*-
#
#  This file is part of BioKit software
#
#  Copyright (c) 2016 - BioKit Development Team
#
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import time
import abc
import select
import sys
import inspect

from io import StringIO
from subprocess import Popen, PIPE

import colorlog
_log = colorlog.getLogger(__name__)

from bioconvert.core.benchmark import Benchmark

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
    def __init__(cls, name, bases, classdict):

        def check_ext(ext, io_name):
            """
            Check if the extension is specified correctly.
            I must be a string or a sequence of string, otherwise raise an error
            it should start with a dot. Otherwise fix extension and inject it in the class

            :param ext: the value of the class attribute (input|output)_ext
            :type ext: a string or a list, tuple or set of strings
            :param str io_name: the type of extension, 'input' or output'
            :raise TypeError:  if ext is neither a string nor a sequence of strings
            """
            if isinstance(ext, str):
                if not ext.startswith('.'):
                    ext = '.' + ext
                setattr(cls, '{}_ext'.format(io_name),  (ext, ))
            elif isinstance(ext, (list, tuple, set)):
                if not all([isinstance(one_ext, str) for one_ext in ext]):
                    raise TypeError("each element of the class attribute '{}.{}_ext' "
                                    "must be a string".format(cls, io_name))
                else:
                    if not all([one_ext.startswith('.') for one_ext in ext]):
                        fixed_ext = []
                        for one_ext in ext:
                            if one_ext.startswith('.'):
                                fixed_ext.append(one_ext)
                            else:
                                fixed_ext.append('.' + one_ext)
                        setattr(cls, '{}_ext'.format(io_name), fixed_ext)
            else:
                import sys
                err = "the class attribute '{}.{}_ext' " \
                      "must be specified in the class or subclasses".format(cls.__name__, io_name)
                _log.warning("skip class '{}': {}".format(cls.__name__, err, file=sys.stderr))
                raise TypeError("the class attribute '{}.{}_ext' must be specified "
                                "in the class or subclasses".format(cls.__name__, io_name))
            return True

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
            if '2' not in name:
                raise TypeError("converter name must follow convention input2output")
            input_fmt, output_fmt = name.upper().split('2', 1)
            input_ext = getattr(cls, 'input_ext')
            if check_ext(input_ext, 'input'):
                output_ext = getattr(cls, 'output_ext')
                check_ext(output_ext, 'output')
            setattr(cls, 'input_fmt', input_fmt)
            setattr(cls, 'output_fmt', output_fmt)
            available_conv_meth = inspect.getmembers(cls, is_conversion_method)

            # do not use strip() but split()
            available_conv_meth = [name[0].split("_method_")[1] for name in available_conv_meth]
            setattr(cls, 'available_methods', available_conv_meth)
            _log.debug("class = {}  available_methods = {}".format(cls.__name__, available_conv_meth))


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

    """specify the extensions of the input file, can be a sequence (must be 
    overridden in subclasses)"""
    input_ext = None
    """specify the extensions of the output file, can be a sequence (must be 
    overridden in subclasses)"""
    output_ext = None
    _default_method = None

    def __init__(self, infile, outfile):
        """

        :param str infile: The path of the input file.
        :param str outfile: The path of The output file
        """
        self.infile = infile
        self.outfile = outfile
        from easydev.multicore import cpu_count
        self.threads = cpu_count()

    def __call__(self, *args, threads=None, **kwargs):
        """

        """
        # If method provided, use it 
        method_name = kwargs.get("method", None)
        if method_name:
            del kwargs["method"]

        # If not, but there is one argument, presumably this is
        # the method
        if method_name is None and len(args) == 1:
            method_name = args[0]
        elif method_name is None and len(args) == 0:
            method_name = self.default

        # If not, we need to check the name
        # "dummy" is a method used to evaluate the cost of the
        # execute() method for the benchmark
        if method_name not in self.available_methods + ['dummy']:
            msg = "Method available are {}".format(self.available_methods)
            _log.error(msg)
            raise ValueError(msg)

        _log.info("Executing {} method".format(method_name))
        # reference to the method requested
        method_reference = getattr(self, "_method_{}".format(method_name))

        # call the method itself
        method_reference(*args, threads=None, **kwargs)

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

    def execute(self, cmd, ignore_errors=False, verbose=False):
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

    def boxplot_benchmark(self, N=5, rerun=True, include_dummy=False):
        """Simple wrapper to call :class:`Benchmark` and plot the results"""
        self._benchmark = Benchmark(self, N=N)
        self._benchmark.include_dummy = include_dummy
        self._benchmark.plot(rerun=rerun)

    def _get_default_method(self):
        if self._default_method is None:
            return self.available_methods[0]
        else:
            return self._default_method
    default = property(_get_default_method)
