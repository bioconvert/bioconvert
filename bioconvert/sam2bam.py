"""Convert :term:`SAM` file to :term:`BAM` file"""
from bioconvert import ConvBase


class SAM2BAM(ConvBase):
    """
    command used::
        samtools view -Sbh
    """
    input_ext = '.sam'
    output_ext ='.bam'

    def _method_samtools(self, *args, **kwargs):
        """
        Do the conversion  sorted :term`SAM` -> :term:'BAM`
        The result of the conversion is stored in the outputfile 
        :return: the standard output
        :rtype: :class:`io.StringIO` object.
        """
        # -S means ignored (input format is auto-detected)
        # -b means output to BAM format
        # -h means include header in SAM output
        cmd = "samtools view -Sbh {} > {}".format(self.infile, self.outfile)
        self.execute(cmd)

