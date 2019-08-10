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
        self._ext_registry = {}
        self._fmt_registry = {}
        self._fill_registry(bioconvert.__path__)
        self._build_path_dict()

    def _fill_registry(self, path, target=None, including_not_available_converter=False):
        """
        Explore the directory converters to discover all converter classes
        (a concrete class which inherits from :class:`ConvBase`)
        and fill the register with the input format and output format
        associated to this converter.

        This is called in the constructor once with
        including_not_available_converter set to False and called at any time 
        to :meth:`get_all_conversions` with including_not_available_converter
        set to True.

        :param str path: the path of a directory to explore (not recursive)
        :param str target:
        :param bool including_not_available_converter:
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
                        format_pair = (converter.input_fmt, converter.output_fmt)
                        target[(format_pair)] = converter
                        # have all the combinaisons between the extensions of 
                        # output formats of the convertes
                        combo_input_ext = tuple(itertools.product(*converter.input_ext))
                        # have all the combinaisons between the extensions of output 
                        # formats of the convertes
                        combo_output_ext = tuple(itertools.product(*converter.output_ext))
                        all_ext_pair = tuple(itertools.product(combo_input_ext,(combo_output_ext)))
                        for ext_pair in all_ext_pair:
                            if len(converter.available_methods) == 0 and not including_not_available_converter:
                                _log.warning("converter '{}' for {} -> {} was not added as no method is available"
                                             .format(converter_name, *ext_pair))
                            else:
                                self.set_ext(ext_pair, converter)

    def _build_path_dict(self):
        """
        Construct dictionaries of dictionaries containing shortest paths
        from one format to another.
        """
        from networkx import DiGraph, all_pairs_shortest_path
        # all_pairs_shortest_path yields pairs (n, d) where
        # * n is a node
        # * d is a dict where
        #     * keys are other nodes
        #     * values are shortest paths (i.e. lists of nodes)
        #       from n to these nodes
        # Converting this to dict results in a dict of dicts where
        # * the first key is the source node
        # * the second key is the destination node
        # * the value is a shortest path between these nodes.

        # self._path_dict[in_fmt][out_fmt] = [in_fmt, ..., out_fmt]
        self._path_dict = dict(all_pairs_shortest_path(
            # Directed graph of available in_fmt -> out_fmt conversions
            DiGraph(self.get_conversions())))

        # self._path_dict_ext[in_ext][out_ext] = [in_ext, ..., out_ext]
        self._path_dict_ext = dict(all_pairs_shortest_path(
            # Directed graph of available in_ext -> out_ext conversions
            DiGraph(self.get_conversions_from_ext())))

    def conversion_path(self, input_fmt, output_fmt):
        """
        Return a list of conversion steps to get from input and
        output formats

        :param tuple input_fmt:
        :param tuple output_fmt:

        Each step in the list is a pair of formats.
        """
        try:
            fmt_steps = self._path_dict[input_fmt][output_fmt]
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
            raise KeyError('an other converter already exists for {} -> {}'
                           .format("_".join(format_pair[0]),"_".join(format_pair[1])))
        self._fmt_registry[format_pair] = convertor

    def _check_input_ext(self, ext_pair):
        assert len(ext_pair) == 2, "parameter must be a tuple with 2 items"
        assert isinstance(ext_pair[0], tuple), "first item must be a tuple"
        assert isinstance(ext_pair[1], tuple), "second item must be a tuple"

    def set_ext(self, ext_pair, convertor):
        """
        Register new convertor from input extension and output extension
        in a list. We can have a list of multiple convertors for one
        ext_pair.

        :param tuple ext_pair: tuple containing the input extensions and the
            output extensions e.g. ( ("fastq",) , ("fasta") )
        :param convertor: the convertor which handle the conversion
                          from input_ext -> output_ext
        :type convertor: list of :class:`ConvBase` object
        """
        self._check_input_ext(ext_pair)
        if ext_pair in self._ext_registry:
            self._ext_registry[ext_pair].append(convertor)
        else:
            self._ext_registry[ext_pair] = [convertor]

    def __getitem__(self, format_pair):
        """
        :param format_pair: the input format, the output format
        :type format_pair: tuple of 2 strings
        :return: an object of subclass o :class:`ConvBase`
        """
        format_pair = (format_pair[0], format_pair[1])
        return self._fmt_registry[format_pair]

    def get_ext(self, ext_pair):
        """
        Copy the registry into a dict that behaves like a list
        to be able to have multiple values for a single key
        and from a key have all converter able to do the conversion
        from the input extension to the output extension.

        :param ext_pair: the input extension, the output extension
        :type ext_pair: tuple of 2 strings
        :return: list of objects of subclass o :class:`ConvBase`
        """
        self._check_input_ext(ext_pair)
        return self._ext_registry[ext_pair]

    def __contains__(self, format_pair):
        """
        Can use membership operation on registry to test if a converter
        to go form input format to output format exists.

        :param format_pair: the input format, the output format
        :type format_pair: tuple (or list) of 2 items. The items must be a
            string or a tuple/list of strings.
        :return: True if format_pair is in registry otherwise False.
        """

        # make sure input is tuple of 2 items
        if isinstance(format_pair, (tuple, list)) is False:
            raise ValueError("input argument must be a tuple or list of 2 items")

        if isinstance(format_pair, list):
            format_pair = tuple(format_pair)

        # make sure we have a pair
        if len(format_pair) != 2:
            raise ValueError( "input must have 2 items")

        # make sure each item is a string or tuple/list and convert into tuples
        # first item
        if isinstance(format_pair[0], str):
            format_pair = ( (format_pair[0], ), format_pair[1])
        elif isinstance(format_pair[0], list):
            format_pair = ( tuple(format_pair[0]), format_pair[1])

        # second item
        if isinstance(format_pair[1], str):
            format_pair = ( format_pair[0], (format_pair[1],))
        elif isinstance(format_pair[1], list):
            format_pair = ( format_pair[0], tuple(format_pair[1]))

        # make sure it is upper case
        for item in format_pair[1]:
            item = item.upper()
        for item in format_pair[0]:
            item = item.upper()

        return format_pair in self._fmt_registry

    def __iter__(self):
        """
        Make registry iterable
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

    def get_converters_names(self):
        """
        :return: a generator that allows to get the name of the converter
                 from the subclass (ConvBase object)
        :rtype: generator

        """
        for converter in self._fmt_registry.values():
            yield converter.__name__.lower()

    def get_conversions_from_ext(self):
        """
        :return: a generator which allow to iterate on all available conversions
                 a conversion is encoded by a tuple of
                 2 strings (input extension, output extension)
        :rtype: generator
        """
        for conv in self._ext_registry:
            yield conv

    def get_all_conversions(self):
        """
        :return: a generator which allow to iterate on all available 
            conversions and their availability;  a conversion is encoded 
            by a tuple of 2 strings (input format, output format)
        :retype: generator (input format, output format, status)
        """
        all_converter = {}
        self._fill_registry(bioconvert.__path__, target=all_converter,
            including_not_available_converter=True)

        for i, o in all_converter:
            yield i, o, (i, o) in self._fmt_registry and len(self._fmt_registry[(i, o)].available_methods) > 0

    def conversion_exists(self, input_fmt, output_fmt, allow_indirect=False):
        """
        :param str input_fmt: the input format
        :param str output_fmt: the output format
        :param boolean allow_indirect: whether to count indirect conversions
        :return: True if a converter which transform input_fmt into output_fmt exists
        :rtype: boolean
        """
        input_fmt = tuple([x.upper() for x in input_fmt])
        output_fmt = tuple([x.upper() for x in output_fmt])

        return ((input_fmt, output_fmt) in self._fmt_registry
                or (allow_indirect
                    and len(self.conversion_path(input_fmt, output_fmt))))

    def get_info(self):
        converters = set([self[this] for this in self._fmt_registry])
        data = {}
        for converter in converters:
            data[converter] = len(converter.available_methods)
        return data

    def iter_converters(self, allow_indirect: bool = False):
        """

        :param bool allow_indirect: also return indirect conversion
        :return: a generator to iterate over (in_fmt, out_fmt, converter 
            class when direct, path when indirect)
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

    def __str__(self):
        data = self.info()
        C = data["converters"]
        F = data["formats"]
        M = data["methods"]
        txt = "Number of formats: {}".format(F)
        txt +="\n" + "Number of converters: {}".format(C)
        txt += "\n" + "Number of methods : {}".format(M)
        return txt

    def info(self):
        info = self.get_info()
        _converters = [x for x in info.items()]
        _data = [info[k] for k,v in info.items()]
        C = len(_converters)
        M = sum(_data)

        F = len(set([x for items in self.get_all_conversions() for x in items]))

        return {
            "formats": F, 
            "converters": C, 
            "methods": M,
            "methods_per_converter": round(float(M)/C,2)}
