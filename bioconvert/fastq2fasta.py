from Bio import SeqIO
from bioconvert import ConvBase


class Fastq2Fasta(ConvBase):
    """
    Convert :term:`FASTQ` to :term:`FASTA`
    """

    input_ext = ['.fastq', '.fq']
    output_ext = '.fasta'


    def __init__(self, inputfile, outputfile):
        """
        :param str infile: The path to the input FASTA file. 
        :param str outfile: The path to the output file
        """
        self.inputfile = inputfile
        self.outputfile = outputfile
        self._default_method = "biopython"

    def _method_biopython(self, *args, **kwargs):
        records = SeqIO.parse(self.inputfile, 'fastq')
        SeqIO.write(records, self.outputfile, 'fasta')

    """def _method_python(self):
        with open(self.infile, "r") as fin:
            with open(self.outfile, "w") as fout:
    """

    def _method_seqtk(self, *args, **kwargs):
        cmd = "seqtk seq -A {} > {}".format(self.inputfile, self.outputfile)
        self.execute(cmd)
