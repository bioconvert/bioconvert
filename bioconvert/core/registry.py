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
"""Main bioconvert registry that fetches automatically the relevant converter"""
import inspect
import itertools
import pkgutil
import importlib
import colorlog

import bioconvert

_log = colorlog.getLogger(__name__)


__all__ = ['Registry']


class Registry(object):
    """class to centralise information about available conversions

    ::

        from bioconvert.core.registry import Registry
        r = Registry()
        r.conversion_exists("BAM", "BED")
        r.info()  # returns number of available methods for each converter

        conv_class = r[(".bam", ".bed")]
        converter = conv_class(input_file, output_file)
        converter.convert()

    """
    def __init__(self):
        # self._ext_registry = {}
        self._fmt_registry = {}
        self._fill_registry(bioconvert.__path__)

    def _fill_registry(self, path):
        """
        Explore the directory converters to discover all converter classes
        (a concrete class which inherits from :class:`ConvBase`)
        and fill the register with the input format and output format associated to this
        converter

        :param str path: the path of a directory to explore (not recursive)
        """
        def is_converter(item):
            """Check if a module is a converter"""
            obj_name, obj = item
            if not inspect.isclass(obj):
                return False

            # Note that on some Python version, the isabstract is buggy.
            # Therefore, the isastract does not return False for ConvBase
            # hence the additional check (obj_name in ["ConvBase"])
            return issubclass(obj, bioconvert.ConvBase) \
                    and not inspect.isabstract(obj) \
                    and obj_name not in ["ConvBase"]

        modules = pkgutil.iter_modules(path=path)
        for _, module_name, *_ in modules:
            if module_name != '__init__':
                try:
                    module = importlib.import_module("bioconvert." + module_name)
                except (ImportError, TypeError) as err:
                    _log.warning("skip module '{}': {}".format(module_name, err))
                    continue

                converters = inspect.getmembers(module)
                converters = [c for c in converters if is_converter(c)]
                for converter_name, converter in converters:
                    if converter is not None:
                        # the registry is no more based on extension but on format
                        # all_conv_path = itertools.product(converter.input_ext, converter.output_ext)
                        # for conv_path in all_conv_path:
                        #     self[conv_path] = converter
                        conv_path = (converter.input_fmt, converter.output_fmt)
                        _log.debug("add converter '%s' for %s -> %s", converter_name, *conv_path)
                        self[conv_path] = converter

    def __setitem__(self, conv_path, convertor):
        """
        Register new convertor from input format to output format.

        :param conv_path: the input format, the output format
        :type conv_path: tuple of 2 strings
        :param convertor: the convertor which handle the conversion from input_fmt -> output_fmt
        :type convertor: :class:`ConvBase` object
        """
        if conv_path in self._fmt_registry:
            raise KeyError('an other converter already exist for {} -> {}'.format(*conv_path))
        self._fmt_registry[conv_path] = convertor

    def __getitem__(self, conv_path):
        """
        :param conv_path: the input format, the output format
        :type conv_path: tuple of 2 strings
        :return: an object of subclass o :class:`ConvBase`
        """
        return self._fmt_registry[conv_path]

    def __contains__(self, conv_path):
        """
        Can use membership operation on registry to test if the it exist a converter to
        go form input format to output format.

        :param conv_path: the input format, the output format
        :type conv_path: tuple of 2 strings
        :return: True if conv_path is in registry otherwise False.
        """
        return conv_path in self._fmt_registry

    def __iter__(self):
        """
        make registry iterable through conv_path (str input format, str output format)
        """
        for path in self._fmt_registry:
            yield path

    # this function is no more needed because __setitem__ takes the place
    # def set_fmt_conv(self, in_fmt, out_fmt, converter):
    #     """
    #     Create an entry in the registry for (in_fmt, out_fmt) and the corresponding converter

    #     :param str in_fmt: the output format
    #     :param str out_fmt: the output format
    #     :param converter: the converter able to convert in_fmt into out_fmt
    #     :type converter:  :class:`BaseConv` concrete class
    #     :return: None
    #     """
    #     self._fmt_registry[(in_fmt, out_fmt)] = converter

    def get_conversions(self):
        """
        :return: a generator which allow to iterate on all available conversions
                 a conversion is encoded by a tuple of 2 strings (input format, output format)
        :retype: generator
        """
        for conv in self._fmt_registry:
            yield conv

    def conversion_exists(self, in_fmt, out_fmt):
        """
        :param str in_fmt: the input format
        :param str out_fmt: the output format
        :return: True if it exists a converter which transform in_fmt into out_fmt
        :rtype: boolean
        """
        in_fmt = in_fmt.upper()
        out_fmt = out_fmt.upper()
        return (in_fmt, out_fmt) in self._fmt_registry

    def get_info(self):
        converters = set([self[this] for this in self._fmt_registry])
        data = {}
        for converter in converters:
            data[converter] = len(converter.available_methods)
        return data

