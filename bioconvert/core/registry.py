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
# import itertools
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
        self._build_path_dict()

    def _fill_registry(self, path, target=None, including_not_available_converter=False):
        """
        Explore the directory converters to discover all converter classes
        (a concrete class which inherits from :class:`ConvBase`)
        and fill the register with the input format and output format
        associated to this converter

        :param str path: the path of a directory to explore (not recursive)
        """

        target = self if target is None else target

        def is_converter(item):
            """Check if a module is a converter"""
            obj_name, obj = item
            if not inspect.isclass(obj):
                return False

            # Note that on some Python version, the isabstract is buggy.
            # Therefore, the isabstract does not return False for ConvBase
            # hence the additional check (obj_name in ["ConvBase"])
            return (issubclass(obj, bioconvert.ConvBase)
                    and not inspect.isabstract(obj)
                    and obj_name not in ["ConvBase"])

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
                        # all_format_pair = itertools.product(
                        #     converter.input_ext, converter.output_ext)
                        # for format_pair in all_format_pair:
                        #     self[format_pair] = converter
                        format_pair = (converter.input_fmt, converter.output_fmt)
                        if len(converter.available_methods) == 0 and not including_not_available_converter:
                            _log.warning("converter '%s' for %s -> %s was not added as no method is available",
                                         converter_name, *format_pair)
                        else:
                            _log.debug("add converter '%s' for %s -> %s",
                                       converter_name, *format_pair)
                            target[format_pair] = converter

    def _build_path_dict(self):
        """
        Constructs a dictionary containing shortest paths
        from one format to another.
        """
        from networkx import DiGraph, all_pairs_shortest_path
        self._path_dict = dict(all_pairs_shortest_path(
            DiGraph(self.get_conversions())))

    # TODO: Should we use a format_pair instead of two strings?
    def conversion_path(self, in_fmt, out_fmt):
        """
        Returns a list of conversion steps to get from *in_fmt* to *out_fmt*.

        Each step in the list is a pair of formats.
        """
        try:
            fmt_steps = self._path_dict[in_fmt][out_fmt]
        except KeyError:
            fmt_steps = []
        return list(zip(fmt_steps, fmt_steps[1:]))

    def __setitem__(self, format_pair, convertor):
        """
        Register new convertor from input format to output format.

        :param format_pair: the input format, the output format
        :type format_pair: tuple of 2 strings
        :param convertor: the convertor which handle the conversion
                          from input_fmt -> output_fmt
        :type convertor: :class:`ConvBase` object
        """
        if format_pair in self._fmt_registry:
            raise KeyError('an other converter already exists for {} -> {}'.format(*format_pair))
        self._fmt_registry[format_pair] = convertor

    def __getitem__(self, format_pair):
        """
        :param format_pair: the input format, the output format
        :type format_pair: tuple of 2 strings
        :return: an object of subclass o :class:`ConvBase`
        """
        return self._fmt_registry[format_pair]

    def __contains__(self, format_pair):
        """
        Can use membership operation on registry to test if a converter
        to go form input format to output format exists.

        :param format_pair: the input format, the output format
        :type format_pair: tuple of 2 strings
        :return: True if format_pair is in registry otherwise False.
        """
        return format_pair in self._fmt_registry

    def __iter__(self):
        """
        make registry iterable
        through format_pair (str input format, str output format)
        """
        for format_pair in self._fmt_registry:
            yield format_pair

    def get_conversions(self):
        """
        :return: a generator which allow to iterate on all available conversions
                 a conversion is encoded by a tuple of
                 2 strings (input format, output format)
        :retype: generator
        """
        for conv in self._fmt_registry:
            yield conv

    def get_all_conversions(self):
        """
        :return: a generator which allow to iterate on all available conversions and their availability
                 a conversion is encoded by a tuple of
                 2 strings (input format, output format)
        :retype: generator (input format, output format, status)
        """
        all_converter = {}
        self._fill_registry(bioconvert.__path__, all_converter, True)
        for i, o in all_converter:
            yield i, o, (i, o) in self._fmt_registry

    def conversion_exists(self, in_fmt, out_fmt, allow_indirect=False):
        """
        :param str in_fmt: the input format
        :param str out_fmt: the output format
        :param boolean allow_indirect: whether to count indirect conversions
        :return: True if a converter which transform in_fmt into out_fmt exists
        :rtype: boolean
        """
        in_fmt = in_fmt.upper()
        out_fmt = out_fmt.upper()
        # return (in_fmt, out_fmt) in self._fmt_registry
        return ((in_fmt, out_fmt) in self._fmt_registry
                or (allow_indirect
                    and len(self.conversion_path(in_fmt, out_fmt))))

    def get_info(self):
        converters = set([self[this] for this in self._fmt_registry])
        data = {}
        for converter in converters:
            data[converter] = len(converter.available_methods)
        return data

    def get_converter(self, in_fmt, out_fmt):
        """

        :param str in_fmt: the format of the input
        :param str out_fmt:  the format of the output
        :return: the converter which convert in_fmt to to_fmt
        :rtype: a :class:`BaseConv` concrete class o
        """
        in_fmt = in_fmt.upper()
        out_fmt = out_fmt.upper()
        return self._fmt_registry((in_fmt, out_fmt))

    def iter_converters(self, allow_indirect:bool=False):
        """

        :param bool allow_indirect: also return indirect conversion
        :return: a generator to iterate over (in_fmt, out_fmt, converter class when direct, path when indirect)
        :rtype: a generator
        """
        # if allow_indirect:
        for start, stops in self._path_dict.items():
            for stop, path in stops.items():
                if len(path) == 1:
                    pass
                elif len(path) == 2:
                    yield start, stop, self._fmt_registry[(start, stop)], None
                elif allow_indirect:
                    yield start, stop, None, path
        #     return
        # for conv, converter in self._fmt_registry.items():
        #     in_fmt, out_fmt = conv
        #     yield in_fmt, out_fmt, converter, None


