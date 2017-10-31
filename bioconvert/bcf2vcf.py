"""Convert :term:`BCF` file to :term:`VCF` file"""
from bioconvert import ConvBase


class BCF2VCF(ConvBase):
    """

    """
    input_ext = [".bcf"]
    output_ext = [".vcf"]

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile:
        :param str outfile:

        command used::

            bcftools view 
        """
        super(BCF2VCF, self).__init__(infile, outfile, *args, **kargs)

    def _method_bcftools(self, *args, **kwargs):

        # -O, --output-type b|u|z|v Output compressed BCF (b), uncompressed BCF
        # (u), compressed VCF (z), uncompressed VCF (v). Use the -Ou option when
        # piping between bcftools subcommands to speed up performance
        cmd = "bcftools view {} -O v -o {}".format(self.infile, self.outfile)
        self.execute(cmd)



