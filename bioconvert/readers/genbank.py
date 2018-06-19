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

class Genbank:
	""" Genbank file format parser.
	This parser is based on the format description presented on the NCIBI website
	https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
	"""
	authorized_keywords = ["LOCUS", "DEFINITION", "ACCESSION", "VERSION", "KEYWORDS", "SOURCE", "REFERENCE", "COMMENT", "FEATURES", "ORIGIN"]

	def __init__(self, filename):
		self.filename = filename
		self.sequence = {}
		self.parsing_function = None

	def read(self):
		with open(self.filename) as reader:
			for line in reader:
				line = " ".join(line.split())
				
				# Empty line
				if len(line) == 0:
					continue

				# Sequence separator found
				if line.startswith("//"):
					if len(self.sequence):
						yield self.sequence
					self.sequence = {}
					continue

				self.parse_routing(line)

		if len(self.sequence) > 0:
			yield self.sequence


	# --------------- Parser router ---------------

	def parse_routing(self, line):
		# Find the right line parser
		if line.startswith("LOCUS"):
			self.parsing_function = self.parse_locus
		elif line.startswith("DEFINITION"):
			self.parsing_function = self.parse_definition
		elif line.startswith("ACCESSION"):
			self.parsing_function = self.parse_accession
		elif line.startswith("VERSION"):
			self.parsing_function = self.parse_version
		elif line.startswith("KEYWORDS"):
			self.parsing_function = self.parse_keywords
		elif line.startswith("SOURCE"):
			self.parsing_function = self.parse_source
		elif line.startswith("REFERENCE"):
			self.parsing_function = self.parse_reference
		elif line.startswith("COMMENT"):
			self.parsing_function = self.parse_comment
		elif line.startswith("FEATURES"):
			self.parsing_function = self.parse_feature
		elif line.startswith("ORIGIN"):
			self.parsing_function = self.parse_origin
		# TODO: Remove this elif when all the parser is over
		elif line[:line.find(" ")] in Genbank.authorized_keywords:
			self.parsing_function = None
			pass
		
		# Parse the line
		if self.parsing_function != None:
			self.parsing_function(line)
		else:
			print("Impossible to parse line\n{}".format(line))


	# --------------- Sub-parser by keyword ---------------

	def parse_locus(self, line):
		""" Parse the LOCUS line.
		The line is splited into 5 values as follow :
		- Sequence id
		- Sequence length
		- Molecule type (can be multiple values. Please read the genbank documentation for values)
		- Genbank division
		- Modification date
		"""
		# split the locus line
		split = line.split()
		# Add false values for if the number of values is incorect
		for _ in range(len(split), 7):
			split += ["None"]
		
		# Parse values of interest
		id, length = split[1:3]
		type = " ".join(split[4:-2])
		division, date = split[-2:]

		# Create dictionary that will contain all the locus values
		self.sequence["LOCUS"] = {"id":id, "length":length, "mol_type":type, "genbank_div":division, "date":date}

	def parse_definition(self, line):
		""" Parse the DEFINITION section.
		This section contains a text without format restriction.
		"""
		# first call
		if not "DEFINITION" in self.sequence:
			# remove the DEFINITION keyword from the text
			self.sequence["DEFINITION"] = line[line.find(" ")+1:]
		else:
			self.sequence["DEFINITION"] = "{} {}".format(self.sequence["DEFINITION"], line)

	def parse_accession(self, line):
		""" Parse the text as a uniq accession number"""
		split = line.split()

		# First line of the accession
		if split[0] == "ACCESSION":
			# Accession id
			self.sequence["ACCESSION"] = {"id": split[1]}

			# Other values from the accession field
			if len(split) > 2:
				self.sequence["ACCESSION"]["other"] = " ".join(split[2:])
		# In case of multi line
		elif "ACCESSION" in self.sequence and "other" in self.sequence["ACCESSION"]:
			self.sequence["ACCESSION"]["other"] = "{} {}".format(self.sequence["ACCESSION"]["other"], " ".join(split))

	def parse_version(self, line):
		""" Parse the text as the version number.
		This value should begin by the accession id a point and then the sequence version number
		Just after the version id the GI identifier can be parsed too
		"""
		split = line.split()

		# First line of the version
		if split[0] == "VERSION":
			# Version id
			self.sequence["VERSION"] = {"id": split[1]}

			# Geninfo Identifier
			if len(split) >2 and split[2].startswith("GI"):
				self.sequence["VERSION"]["GI"] = split[2][split[2].find(":")+1:]
				del split[2]

			# Other values from the Version field
			if len(split) > 2:
				self.sequence["VERSION"]["other"] = " ".join(split[2:])
		# In case of multi line
		elif "VERSION" in self.sequence and "other" in self.sequence["VERSION"]:
			self.sequence["VERSION"]["other"] = "{} {}".format(self.sequence["VERSION"]["other"], " ".join(split))

	def parse_keywords(self, line):
		""" This field is deprecated
		See the documentation for more info
		"""
		if not line.startswith("KEYWORDS ."):
			self.sequence["KEYWORDS"] = line[line.find(" ")+1:]

	def parse_source(self, line):
		""" Parse the SOURCE label.
		Just after the SOURCE keyword there is a short description.
		The complete description is expected after the optional keyword ORGANISM
		"""
		if line.startswith("SOURCE"):
			self.sequence["SOURCE"] = {"short": line[line.find(" ")+1:]}
		elif line.startswith("ORGANISM"):
			self.parsing_function = self.parse_organism
			self.parse_organism(line)
		else:
			self.sequence["SOURCE"]["short"] = "{} {}".format(self.sequence["SOURCE"]["short"], line)

	def parse_organism(self, line):
		""" Parse the source organism.
		After the ORGANISM keyword the name of the organism is expected.
		The folowing lines are dedicated to the lineage description
		"""
		if line.startswith("ORGANISM"):
			self.sequence["SOURCE"]["ORGANISM"] = {"name": line[line.find(" ")+1:]}
		else:
			if "lineage" in self.sequence["SOURCE"]["ORGANISM"]:
				self.sequence["SOURCE"]["ORGANISM"]["lineage"] = "{} {}".format(self.sequence["SOURCE"]["ORGANISM"]["lineage"], line)
			else:
				self.sequence["SOURCE"]["ORGANISM"]["lineage"] = line

	reference_keywords = ["AUTHORS", "TITLE", "JOURNAL", "PUBMED", "REMARK"]
	def parse_reference(self, line):
		""" Parse the publication sections
		The keyword REFERENCE can appear as much as you want
		"""
		# Add the publication list
		if not "REFERENCE" in self.sequence:
			self.sequence["REFERENCE"] = []

		# create new entry
		if line.startswith("REFERENCE"):
			self.ref = {}
			self.sequence["REFERENCE"].append(self.ref)
			return

		# dispatch new values
		keyword = line[:line.find(" ")]
		if keyword == "AUTHORS":
			self.parsing_function = self.parse_authors
		elif keyword == "TITLE":
			self.parsing_function = self.parse_title
		elif keyword == "JOURNAL":
			self.parsing_function = self.parse_journal
		elif keyword == "PUBMED":
			self.parsing_function = self.parse_pubmed
		elif keyword == "REMARK":
			self.parsing_function = self.parse_remark
		else:
			print("Genbank: impossible to parse the following line:\n{}".format(line))
			return

		# Parse keywords for reference
		self.parsing_function(line)

	def is_in_the_correct_subtree(self, line, expected_keyword):
		""" Verify if is in the right subtree of the reference parser.
		If not, parse with the right subtree and return False
		"""
		keyword = line[:line.find(" ")]
		# if must go in another branch of the reference parsing
		if keyword == expected_keyword:
			return True

		if keyword in Genbank.reference_keywords:
			self.parsing_function = self.parse_reference
			self.parsing_function(line)
			return False

	def parse_authors(self, line):
		if not self.is_in_the_correct_subtree(line, "AUTHORS"):
			return

		# Parse the authors
		if line.startswith("AUTHORS"):
			self.ref["AUTHORS"] = line[line.find(" ")+1:]
		else:
			self.ref["AUTHORS"] = "{} {}".format(self.ref["AUTHORS"], line)


	def parse_title(self, line):
		if not self.is_in_the_correct_subtree(line, "TITLE"):
			return

		# Parse the title
		if line.startswith("TITLE"):
			self.ref["TITLE"] = line[line.find(" ")+1:]
		else:
			self.ref["TITLE"] = "{} {}".format(self.ref["TITLE"], line)

	def parse_journal(self, line):
		if not self.is_in_the_correct_subtree(line, "JOURNAL"):
			return

		# Parse the journal
		if line.startswith("JOURNAL"):
			self.ref["JOURNAL"] = line[line.find(" ")+1:]
		else:
			self.ref["JOURNAL"] = "{} {}".format(self.ref["JOURNAL"], line)

	def parse_pubmed(self, line):
		if not self.is_in_the_correct_subtree(line, "PUBMED"):
			return

		# Parse the pubmed id
		if line.startswith("PUBMED"):
			self.ref["PUBMED"] = line[line.find(" ")+1:]
		else:
			print("Genbank parser: unexpected multiline pubmed id")

	def parse_remark(self, line):
		print(line)
		if not self.is_in_the_correct_subtree(line, "REMARK"):
			return

		# Parse the remark line
		if line.startswith("REMARK"):
			self.ref["REMARK"] = line[line.find(" ")+1:]
		else:
			self.ref["REMARK"] = "{} {}".format(self.ref["REMARK"], line)

	def parse_comment(self, line):
		""" Parse the COMMENT field as a text
		"""
		if line.startswith("COMMENT"):
			self.sequence["COMMENT"] = line[line.find(" ")+1:]
		else:
			self.sequence["COMMENT"] = "{} {}".format(self.sequence["COMMENT"], line)

	def parse_feature(self, line):
		print("Genbank FEATURES parsing not yet implemented")
		pass

	def parse_origin(self, line):
		""" Parse the origin label and next lines to construct the sequence """
		# Define the sequence on the label ORIGIN
		if line.startswith("ORIGIN"):
			self.sequence["ORIGIN"] = ""

		# parse the sequence
		else:
			# Extract the sub sequence from the line
			sub_seq = "".join(line.split()[1:])
			# Add the subsequence to the complete sequence
			self.sequence["ORIGIN"] = "{}{}".format(self.sequence["ORIGIN"], sub_seq)



gbk = Genbank("data/test_genbank.gbk")
nb_seq = 0
for seq in gbk.read():
	nb_seq += 1
	print(seq)

print(nb_seq)
