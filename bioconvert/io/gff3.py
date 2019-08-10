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
import re


class Gff3():
    """Read a gff v3 file

    See the format description at
    https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md


    """
    def __init__(self, filename):
        self.filename = filename

    def read(self):
        """ Read annotations one by one creating a generator """
        with open(self.filename) as reader:
            line = None

            for line in reader:
                # Skip metadata and comments
                if line.startswith("#"):
                    continue

                # Format checking
                split = line.rstrip().split("\t")
                if len(split) < 9:
                    # Wrong line format
                    if len(split) > 0:
                        print("Impossible to read the following line regarding the gff3 specifications")
                        print(line)
                    continue

                annotation = self.process_main_fields(split[0:8])
                annotation["attributes"] = self.process_attributes(split[8])

                yield annotation


    def process_main_fields(self, fields):
        annotation = {}

        # Unique id of the sequence
        annotation["seqid"] = self.decode_small(fields[0])
        # Optional source
        if fields[1] != ".":
            annotation["source"] = self.decode_small(fields[1])
        # Annotation type
        annotation["type"] = self.decode_small(fields[2])
        # Start and stop
        annotation["start"] = int(fields[3])
        annotation["stop"] = int(fields[4])
        # Optional score field
        if fields[5] != ".":
            annotation["score"] = float(fields[5])
        # Strand
        if fields[6] == "+" or fields[6] == "-" or fields[6] == "?" or fields[6] == ".":
            annotation["strand"] = fields[6]
        # Phase
        if fields[7] != ".":
            annotation["phase"] = int(fields[7]) % 3

        return annotation


    def process_attributes(self, text):
        attributes = {}

        # split into mutliple attributes
        split = text.split(";")
        for attr in split:
            #find the separator
            idx = attr.find("=")

            # parse tags and associated values
            value = self.decode_complete(attr[idx+1:])
            if len(value) == 1:
                value = value[0]
            attributes[self.decode_complete(attr[:idx])] = value

        return attributes


    @staticmethod
    def decode_small(text):
        # Tabulation
        text = re.sub("%09", "\t", text)
        # newline
        text = re.sub("%0A", "\n", text)
        # return
        text = re.sub("%0D", "\r", text)
        # percent
        text = re.sub("%25", "%", text)

        return text


    @staticmethod
    def decode_complete(text):
        text = Gff3.decode_small(text)
        # semicolon
        text = re.sub("%3B", ";", text)
        # equals
        text = re.sub("%3D", "=", text)
        # ampersand
        text = re.sub("%26", "&", text)
        # comma
        text = re.sub("%2C", ",", text)

        return text


# if __name__ == "__main__":
#     file = "/home/yoann/Projects/Hub/bioconvert/bioconvert/data/gff3_example.gff"
#     gff = Gff3(file)
#     for annotation in gff.read():
#         print(annotation)
