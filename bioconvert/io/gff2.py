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


class Gff2():
	"""Read a gff v3 file
	See the format description at https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
	"""
	def __init__(self, filename):
		self.filename = filename
		

	def read(self):
		""" Read annotations one by one creating a generator """
		
		with open(self.filename) as reader:
			line = None

			for line in reader:
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
		annotation["seqid"] = fields[0]
		# Optional source
		if fields[1] != ".":
			annotation["source"] = fields[1]
		# Annotation type
		annotation["type"] = fields[2]
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
		split = text.rstrip().split(";")
		for attr in split:
			if len(attr) == 0 or attr == ".":
				continue

			# Remove parasite spaces
			attr = attr.strip()

			# Find the key and value separator
			idx = attr.find(" ")
			# split key from value
			key, value = attr[:idx].strip(), attr[idx+1:].strip()
			# format the value
			if value[0] == "\"":
				value = value[1:-1]
			
			if not key in attributes:
				attributes[key] = value
			else:
				if not isinstance(attributes[key], list):
					attributes[key] = [attributes[key]]
				attributes[key].append(value)

		return attributes

			

# if __name__ == "__main__":
# 	file = "/home/yoann/Projects/Hub/bioconvert/bioconvert/data/gff2_example.gff"
# 	gff = Gff2(file)
# 	for annotation in gff.read():
# 		print(annotation)