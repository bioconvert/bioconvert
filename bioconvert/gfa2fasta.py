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
"""Convert :term:`GFA` format to :term:`FASTA` formats"""
from bioconvert import ConvBase
import colorlog
logger = colorlog.getLogger(__name__)



__all__ = ["GFA2FASTA"]


class GFA2FASTA(ConvBase):
    """Convert sorted :term:`GFA` file into :term:`FASTA` file 

    Available methods: awk, python

    .. plot::

         from bioconvert.gfa2fasta import GFA2FASTA
         from bioconvert import bioconvert_data
         from easydev import TempFile

         with TempFile(suffix=".fasta") as fh:
             infile = bioconvert_data("test_gfa2fasta.gfa")
             convert = GFA2FASTA(infile, fh.name)
             convert.boxplot_benchmark()

    :reference: https://github.com/GFA-spec/GFA-spec/blob/master/GFA-spec.md

    .. seealso:: bioconvert.simulator.gfa

    """
    input_ext = ['.gfa']
    output_ext = ['.fasta', ".fa"]

    def __init__(self, infile, outfile):
        """
        :param str infile: The path to the input BAM file. **It must be sorted**.
        :param str outfile: The path to the output file
        """
        super().__init__(infile, outfile)
        self._default_method = "python"

    def _method_awk(self, *args, **kwargs):
        """

        :return: the standard output
        :rtype: :class:`io.StringIO` object.

        .. note:: this method fold the sequence to 80 characters
        """
        # Note1: since we use .format, we need to escape the { and } characters
        # Note2: the \n need to be escaped for Popen to work
        cmd = """awk '/^S/{{print ">"$2"\\n"$3}}' {} | fold > {}""".format(self.infile, self.outfile)
        self.execute(cmd)

    def _method_python(self, *args, **kwargs):
        with open(self.infile, "r") as fin:
            with open(self.outfile, "w") as fout:
                for line in fin.readlines():
                    if line.startswith("S"):
                        field, name, sequence = line.split()
                        fout.write(">{}\n{}\n".format(name, sequence))

