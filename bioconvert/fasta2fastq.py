"""Convert :term:`Fastq` format to :term:`Fastq` formats"""
from bioconvert import ConvBase
from pysam import FastxFile

__all__ = ["Fasta2Fastq"]


class Fasta2Fastq(ConvBase):
    """

    """
    input_ext = ['.fa', '.fas', '.fasta']
    output_ext = ['.fastq', 'fq']

    def __init__(self, infile, outfile):
        """
        :param str infile: The path to the input FASTA file
        :param str outfile: The path to the output FASTQ file
        """
        super().__init__(infile, outfile)
        self._default_method = "v1"

    def _method_v1(self, *args, **kwargs):

        with open(self.outfile, 'w') as fastq_out:
        
            for seq in FastxFile(self.infile):
                fastq_out.write("@{0} {1}\n{2}\n+\n{3}\n".format(seq.name,
                                                                 seq.comment,
                                                                 seq.sequence,
                                                                 len(seq.sequence) * "I"))

                
