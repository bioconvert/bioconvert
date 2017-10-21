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

        self.available_methods = ["samtools", "bedtools"]
        self.default = "samtools"

    def __call__(self, *args, **kwargs):
        """
        do the conversion  sorted :term`BAM` -> :term:'BED`

        :return: the standard output
        :rtype: :class:`io.StringIO` object.
        """
        method = kwargs.get("method", self.default)
        assert method in self.available_methods
        if method == "samtools":
            self._method_samtools(*args, **kwargs)
        elif method == "bedtools":
            self._method_bedtools(*args, **kwargs)

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
