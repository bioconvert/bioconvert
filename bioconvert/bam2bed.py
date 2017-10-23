"""Convert :term:`BAM` format to :term:`BED` formats"""
from bioconvert import ConvBase

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
        self.default = "samtools"

    def __call__(self, *args, **kwargs):
        """
        do the conversion  sorted :term`BAM` -> :term:'BED`

        :return: the standard output
        :rtype: :class:`io.StringIO` object.
        """
        method_name = kwargs.get("method", self.default)
        assert method_name in self.available_methods
        method_reference = getattr(self, "_method_{}".format(method_name))
        method_reference(*args, **kwargs)

    def _method_samtools(self, *args, **kwargs):
        """
        do the conversion  sorted :term`BAM` -> :term:'BED` using samtools

        :return: the standard output
        :rtype: :class:`io.StringIO` object.
        """
        cmd = "samtools depth -aa {} > {}".format(self.infile, self.outfile)
        self.execute(cmd)

    def _method_bedtools(self, *args, **kwargs):
        """
        do the conversion  sorted :term`BAM` -> :term:'BED` using samtools

        :return: the standard output
        :rtype: :class:`io.StringIO` object.
        """
        cmd = "bedtools genomecov -d -ibam {} > {}".format(self.infile, self.outfile)
        self.execute(cmd)
