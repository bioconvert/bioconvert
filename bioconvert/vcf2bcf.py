"""Convert :term:`VCF` file to :term:`BCF` file"""
from biokit.converters.convbase import ConvBase


class VCF2BCF(ConvBase):
    """

    """
    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile:
        :param str outfile:

        command used::

            bcftools view -Sb
        """
        super(VCF2BCF, self).__init__(infile, outfile, *args, **kargs)

    def convert(self):
        # -S means ignored (input format is VCF)
        # -b output BCF instead of VCF
        cmd = "bcftools view -Sb {} >  {}".format(self.infile, self.outfile)
        self.execute(cmd)



