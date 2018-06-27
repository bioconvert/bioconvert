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

class Fasta():
	"""Read a fasta file"""
	def __init__(self, filename):
		self.filename = filename

	def read(self):
		""" Read fasta sequence by sequence creating a generator """
		sequence = {"id":"", "comment":"", "value":""}
		
		with open(self.filename) as reader:
			line = None

			for line in reader:
				# Remove line return and useless spaces
				line = " ".join(line.split())

				# Header case
				if line[0] == ">":
					if len(sequence["id"]) > 0:
						yield sequence

					# create new sequence
					sequence = {"id":"", "comment":"", "value":""}
					# Parse header to split into id and comment
					idx = line.find(" ")
					if idx == -1:
						sequence["id"] = line[1:]
					else:
						sequence["id"] = line[1:idx]
						sequence["comment"] = line[idx+1:]

				# Sequence case
				else:
					sequence["value"] = "{}{}".format(sequence["value"], line)
		
		if len(sequence["id"]) > 0:
			yield sequence
			