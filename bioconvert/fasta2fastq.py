"""Convert :term:`Fastq` format to :term:`Fastq` formats"""
from bioconvert import ConvBase
import colorlog

from bioconvert.core.base import ConvArg
from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)


__all__ = ["Fasta2Fastq"]


class Fasta2Fastq(ConvBase):
    """

    """
    input_ext = ['.fa', '.fas', '.fasta']
    output_ext = ['.fastq', 'fq']
    _default_method = "pysam"

    def __init__(self, infile, outfile):
        """
        :param str infile: The path to the input FASTA file
        :param str outfile: The path to the output FASTQ file
        """
        super().__init__(infile, outfile)

    @requires(python_library="pysam")
    def _method_pysam(self, quality_file=None, *args, **kwargs):
        from pysam import FastxFile
        if quality_file is None:
            _log.warning("No quality file provided. Please use --quality-file")
            with open(self.outfile, 'w') as fastq_out:
                for seq in FastxFile(self.infile):
                    fastq_out.write("@{0} {1}\n{2}\n+\n{3}\n".format(seq.name,
                                                                 seq.comment,
                                                                 seq.sequence,
                                                                 len(seq.sequence) * "I"))
        else: # length must be equal and identifiers sorted similarly
            with open(self.outfile, "w") as fastq_out:
                for seq, qual in zip(FastxFile(self.infile), FastxFile(quality_file)):
                    assert seq.name == qual.name
                    fastq_out.write("@{0} {1}\n{2}\n+\n{3}\n".format(seq.name,
                                                                 seq.comment,
                                                                 seq.sequence,
                                                                 qual.sequence))


    @classmethod
    def get_additional_arguments(cls):
        yield ConvArg(
            names="--quality-file",
            nargs="?",
            default=None,
            help="The path to the quality file.",
        )
