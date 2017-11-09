from Bio import SeqIO
from Bio.SeqIO import FastaIO
from bioconvert import ConvBase
from gatb import Bank


class Fastq2Fasta(ConvBase):
    """
    Convert :term:`FASTQ` to :term:`FASTA`
    """

    input_ext = ['.fastq', '.fq']
    output_ext = '.fasta'

    @staticmethod
    def unwrap_fasta(infile, outfile):
        """
        This method reads fasta sequences from *infile*
        and writes them unwrapped in *outfile*.
        :param str infile: The path to the input FASTA file.
        :param str outfile: The path to the output file.
        """
        with open(outfile, "w") as fasta_out:
            FastaIO.FastaWriter(fasta_out, wrap=None).write_file(
                SeqIO.parse(infile, 'fasta'))

    def __init__(self, infile, outfile):
        """
        :param str infile: The path to the input FASTA file.
        :param str outfile: The path to the output file.
        """
        self.inputfile = infile
        self.outputfile = outfile
        self._default_method = "biopython"

    def _method_biopython(self, *args, **kwargs):
        records = SeqIO.parse(self.inputfile, 'fastq')
        SeqIO.write(records, self.outputfile, 'fasta')

    def _method_seqtk(self, *args, **kwargs):
        cmd = "seqtk seq -A {} > {}".format(self.inputfile, self.outputfile)
        self.execute(cmd)

    def _method_GATB(self, *args, **kwargs):
        with open(self.outputfile, "w") as fasta:
            for record in Bank(self.inputfile):
                fasta.write(">{}\n{}\n".format(
                    record.comment.decode("utf-8"),
                    record.sequence.decode("utf-8")))
