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
        super().__init__(infile, outfile)
        self._default_method = "biopython"

    def _method_biopython(self, *args, **kwargs):
        records = SeqIO.parse(self.infile, 'fastq')
        SeqIO.write(records, self.outfile, 'fasta')

    def _method_seqtk(self, *args, **kwargs):
        cmd = "seqtk seq -A {} > {}".format(self.infile, self.outfile)
        self.execute(cmd)

    def _method_GATB(self, *args, **kwargs):
        with open(self.outfile, "w") as fasta:
            for record in Bank(self.infile):
                fasta.write(">{}\n{}\n".format(
                    record.comment.decode("utf-8"),
                    record.sequence.decode("utf-8")))

    def _method_awk(self, *args, **kwargs):
        # Note1: since we use .format, we need to escape the { and } characters
        # Note2: the \n need to be escaped for Popen to work
        awkcmd = """awk '{{if(NR%4==1) {{printf(">%s\\n",substr($0,2));}} else if(NR%4==2) print;}}' """
        cmd = "{} {} > {}".format(awkcmd, self.infile, self.outfile)
        self.execute(cmd)

    def _method_mawk(self, *args, **kwargs):
        """This variant of the awk method uses mawk, a lighter and faster implementation of awk."""
        # Note1: since we use .format, we need to escape the { and } characters
        # Note2: the \n need to be escaped for Popen to work
        awkcmd = """mawk '{{if(NR%4==1) {{printf(">%s\\n",substr($0,2));}} else if(NR%4==2) print;}}' """
        cmd = "{} {} > {}".format(awkcmd, self.infile, self.outfile)
        self.execute(cmd)

    def _method_bioawk(self, *args, **kwargs):
        """This method uses bioawk, an implementation with convenient bioinformatics parsing features."""
        awkcmd = """bioawk -c fastx '{{print ">"$name" "$comment"\\n"$seq}}'"""
        cmd = "{} {} > {}".format(awkcmd, self.infile, self.outfile)
        self.execute(cmd)
