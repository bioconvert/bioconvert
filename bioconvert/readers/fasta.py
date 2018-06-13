

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
				# Remove line return
				if line[-1] == "\n":
					line = line[:-1]

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