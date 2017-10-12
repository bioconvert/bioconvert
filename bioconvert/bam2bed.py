"""Convert :term:`BAM` format to :term:`BED` formats"""
from .base import ConvBase

__all__ = ["Bam2Bed"]


class Bam2Bed(ConvBase):
    """
    Convert sorted :term:`BAM` file into :term:`BED` file ::

        samtools depth -aa INPUT > OUTPUT

    """
    input_ext = ['.bam']
    output_ext = '.bed'

    def __init__(self, infile, outfile):
        """
        :param str infile: The path to the input BAM file. **It must be sorted**.
        :param str outfile: The path to the output file
        """
        super().__init__(infile, outfile)

    def __call__(self, *args, **kwargs):
        """
        do the conversion  sorted :term`BAM` -> :term:'BED`

        :return: the standard output
        :rtype: :class:`io.StringIO` object.
        """
        cmd = "samtools depth -aa {} > {}".format(self.infile, self.outfile)
        self.execute(cmd)

