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

		return True

	def parse_authors(self, line):
		if not self.is_in_the_correct_subtree(line, "AUTHORS"):
			return

		# Parse the authors
		if line.startswith("AUTHORS"):
			self.ref["AUTHORS"] = line[line.find(" ")+1:]
		else:
			self.ref["AUTHORS"] = "{} {}".format(self.ref["AUTHORS"], line)


	def parse_title(self, line):
		if not self.is_in_the_correct_subtree(line, "TITLE"):
			return

		# Parse the title
		if line.startswith("TITLE"):
			self.ref["TITLE"] = line[line.find(" ")+1:]
		else:
			self.ref["TITLE"] = "{} {}".format(self.ref["TITLE"], line)

	def parse_journal(self, line):
		if not self.is_in_the_correct_subtree(line, "JOURNAL"):
			return

		# Parse the journal
		if line.startswith("JOURNAL"):
			self.ref["JOURNAL"] = line[line.find(" ")+1:]
		else:
			self.ref["JOURNAL"] = "{} {}".format(self.ref["JOURNAL"], line)

	def parse_pubmed(self, line):
		if not self.is_in_the_correct_subtree(line, "PUBMED"):
			return

		# Parse the pubmed id
		if line.startswith("PUBMED"):
			self.ref["PUBMED"] = line[line.find(" ")+1:]
		else:
			print("Genbank parser: unexpected multiline pubmed id")

	def parse_remark(self, line):
		if not self.is_in_the_correct_subtree(line, "REMARK"):
			return

		# Parse the remark line
		if line.startswith("REMARK"):
			self.ref["REMARK"] = line[line.find(" ")+1:]
		else:
			self.ref["REMARK"] = "{} {}".format(self.ref["REMARK"], line)

	def parse_comment(self, line):
		""" Parse the COMMENT field as a text
		"""
		if line.startswith("COMMENT"):
			self.sequence["COMMENT"] = line[line.find(" ")+1:]
		else:
			self.sequence["COMMENT"] = "{} {}".format(self.sequence["COMMENT"], line)

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

	# --------------- Feature parsing ---------------

	# List of feature keywords
	features = ["assembly_gap", "C_region", "CDS", "centromere", "D-loop", "D_segment", "exon", "gap", "gene", "iDNA", "intron", "J_segment", "mat_peptide", "misc_binding", "misc_difference", "misc_feature", "misc_recomb", "misc_RNA", "misc_structure", "mobile_element", "modified_base", "mRNA", "ncRNA", "N_region", "old_sequence", "operon", "oriT", "polyA_site", "precursor_RNA", "prim_transcript", "primer_bind", "propeptide", "protein_bind", "regulatory", "repeat_region", "rep_origin", "rRNA", "S_region", "sig_peptide", "source", "stem_loop", "STS", "telomere", "tmRNA", "transit_peptide", "tRNA", "unsure", "V_region", "V_segment", "variation", "3'UTR", "5'UTR"]
	# List of qualifiers
	qualifiers = ["allele", "altitude", "anticodon", "artificial_location", "bio_material", "bound_moiety", "cell_line", "cell_type", "chromosome", "citation", "clone", "clone_lib", "codon_start", "collected_by ", "collection_date ", "compare", "country", "cultivar", "culture_collection", "db_xref", "dev_stage", "direction", "EC_number", "ecotype", "environmental_sample", "estimated_length", "exception", "experiment", "focus", "frequency", "function", "gap_type", "gene", "gene_synonym", "germline", "haplogroup", "haplotype", "host", "identified_by ", "inference", "isolate", "isolation_source", "lab_host", "lat_lon ", "linkage_evidence", "locus_tag", "macronuclear", "map", "mating_type", "mobile_element_type", "mod_base", "mol_type", "ncRNA_class", "note", "number", "old_locus_tag", "operon", "organelle ", "organism", "partial", "PCR_conditions", "PCR_primers", "phenotype", "plasmid", "pop_variant", "product", "protein_id", "proviral", "pseudo", "pseudogene", "rearranged", "recombination_class", "regulatory_class", "replace", "ribosomal_slippage", "rpt_family", "rpt_type", "rpt_unit_range", "rpt_unit_seq", "satellite", "segment", "serotype", "serovar", "sex", "specimen_voucher", "standard_name", "strain", "sub_clone", "submitter_seqid", "sub_species", "sub_strain", "tag_peptide", "tissue_lib", "tissue_type", "transgenic", "translation", "transl_except", "transl_table", "trans_splicing ", "type_material", "variety"]

	def parse_feature(self, line):
		""" Parse the features.
		The list of feature used is extracted from the appendix II of the oficial ref.
		For more informations: http://www.insdc.org/documents/feature_table.html#7.2
		The list of qualifier are also extracted from the official documetation, appendix III
		More informations: http://www.insdc.org/documents/feature_table.html#7.3
		"""
		keyword = line[:line.find(" ")]

		# Global features keyword
		if keyword == "FEATURES":
			self.sequence["FEATURES"] = []

		# Feature keyword
		elif keyword in Genbank.features:
			self.current_feature = {"type":keyword, "position":line[line.find(" ")+1:]}
			self.sequence["FEATURES"].append(self.current_feature)

		# Qualifier
		elif line[1: line.find("=") if "=" in line else len(line)] in Genbank.qualifiers:
			# init values
			qualifier = line[1:]
			value = None

			# parse qualifier and value if the operator is dual
			is_dual_operator = "=" in keyword

			if is_dual_operator:
				qualifier, value = qualifier.split("=")
				# remove quotes from strings
				if value[0] == value[-1] == "\"":
					value = value[1:-1]

			# First time that the qualifier is encountered
			if not qualifier in self.current_feature:
				self.current_feature[qualifier] = value
			else:
				# If there is more than 1 element for this qualifier
				if isinstance(self.current_feature[qualifier], list):
					self.current_feature[qualifier].append(value)
				# Transform a value into a vector of values
				else:
					self.current_feature[qualifier] = [self.current_feature[qualifier], value]

		# End of the FEATURES part
		elif keyword in Genbank.authorized_keywords:
			self.parsing_function = self.parse_routing
			self.parsing_function(line)

		# Unknown line
		else:
			print("Genbank: impossible to parse the following line:")
			print(line)



# gbk = Genbank("data/test_genbank.gbk")
# nb_seq = 0
# for seq in gbk.read():
# 	nb_seq += 1
# 	print(seq)

# print(nb_seq)
