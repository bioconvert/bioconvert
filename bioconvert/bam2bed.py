"""Convert :term:`BAM` format to :term:`BED` formats"""
from bioconvert import ConvBase


__all__ = ["BAM2BED"]


class BAM2BED(ConvBase):
    """Convert sorted :term:`BAM` file into :term:`BED` file 

    Available methods:

    - samtools::

        samtools depth -aa INPUT > OUTPUT

    - bedtools::

        bedtools genomecov -d -ibam INPUT > OUTPUT


    .. plot::

         from bioconvert.bam2bed import BAM2BED
         from bioconvert import bioconvert_data
         from easydev import TempFile

         with TempFile(suffix=".bed") as fh:
             infile = bioconvert_data("test_measles.sorted.bam")
             convert = BAM2BED(infile, fh.name)
             convert.boxplot_benchmark()

    """
    input_ext = ['.bam']
    output_ext = '.bed'

    def __init__(self, infile, outfile):
        """
        :param str infile: The path to the input BAM file. **It must be sorted**.
        :param str outfile: The path to the output file
        """
        super().__init__(infile, outfile)
        self._default_method = "samtools"

    def _method_samtools(self, *args, **kwargs):
        """
        do the conversion sorted :term`BAM` -> :term:'BED` using samtools

        :return: the standard output
        :rtype: :class:`io.StringIO` object.
        """
        cmd = "samtools depth -aa {} > {}".format(self.infile, self.outfile)
        self.execute(cmd)

    def _method_bedtools(self, *args, **kwargs):
        """
        do the conversion sorted :term`BAM` -> :term:'BED` using bedtools

        :return: the standard output
        :rtype: :class:`io.StringIO` object.
        """
        cmd = "bedtools genomecov -d -ibam {} > {}".format(self.infile, self.outfile)
        self.execute(cmd)
