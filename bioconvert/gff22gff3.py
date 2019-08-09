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

"""Convert :term:`GFF2` to :term:`GFF3` format"""
from bioconvert import ConvBase
from bioconvert.core.decorators import requires_nothing
from bioconvert import compressor
from bioconvert.io.gff2 import Gff2


__all__ = ["GFF22GFF3"]


class GFF22GFF3(ConvBase):
    """Convert :term:`GFF2` to :term:`GFF3`"""

    _default_method = "bioconvert"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: 
        :param str outfile:

        Method available is pure Python.

        """
        super(GFF22GFF3, self).__init__(infile, outfile, *args, **kargs)

    @requires_nothing
    @compressor
    def _method_bioconvert(self, *args, **kwargs):
        """ This method is a basic mapping of the 9th column of gff2 to gff3.
        Other methods with smart translations must be created for specific usages.
        There is no good solution for this translation.
        """
        gff2_reader = Gff2(self.infile)

        with open(self.outfile, "w") as gff3_writer:
            for annotation in gff2_reader.read():
                # Write the 8 first columns
                gff3_writer.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(
                    annotation["seqid"],
                    annotation["source"] if "source" in annotation else ".", # optional sourse
                    annotation["type"], # should be verified for ontology matching
                    annotation["start"],
                    annotation["stop"],
                    annotation["score"] if "score" in annotation else ".", # Score
                    annotation["strand"] if "strand" in annotation else ".", # starnd
                    annotation["phase"] if "phase" in annotation else "." # phase
                ))
                # Write the 9th column using new standards but without smart translations
                if len(annotation["attributes"]) > 0:
                    attributes = []
                    for key, value in annotation["attributes"].items():
                        # join multiple values
                        if isinstance(value, list):
                            value = ",".join(value)
                        # add the value
                        attributes.append("{}={}".format(key, value))

                    # Write the file
                    gff3_writer.write(";".join(attributes) + "\n")
                else:
                    gff3_writer.write(".\n")
