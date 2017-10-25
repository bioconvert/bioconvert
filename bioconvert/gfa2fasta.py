"""Convert :term:`BAM` format to :term:`BED` formats"""
from bioconvert import ConvBase


__all__ = ["GFA2FASTA"]


class GFA2FASTA(ConvBase):
    """Convert sorted :term:`GFA` file into :term:`FASTA` file 

    Available methods:

    .. plot::

         from bioconvert.gfa2fasta import GFA2FASTA
         from bioconvert import bioconvert_data
         from easydev import TempFile

         with TempFile(suffix=".fasta") as fh:
             infile = bioconvert_data("test_v1.gfa")
             convert = GFA2FASTA(infile, fh.name)
             convert.boxplot_benchmark()

    :reference: https://github.com/GFA-spec/GFA-spec/blob/master/GFA-spec.md
    """
    input_ext = ['.gfa']
    output_ext = ['.fasta', ".fa"]

    def __init__(self, infile, outfile):
        """
        :param str infile: The path to the input BAM file. **It must be sorted**.
        :param str outfile: The path to the output file
        """
        super().__init__(infile, outfile)
        self._default_method = "awk"

    def _method_awk(self, *args, **kwargs):
        """
        do the conversion  sorted :term`BAM` -> :term:'BED` using samtools

        :return: the standard output
        :rtype: :class:`io.StringIO` object.
        """
        # Note1: since we use .format, we need to escape the { and } characters
        # Note2: the \n need to be escaped for Popen to work
        cmd = """awk '/^S/{{print ">"$2"\\n"$3}}' {} | fold > {}""".format(self.infile, self.outfile)
        self.execute(cmd)

