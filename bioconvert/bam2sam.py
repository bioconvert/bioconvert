"""Convert :term:`SAM` file to :term:`BAM` file"""
from bioconvert import ConvBase


class BAM2SAM(ConvBase):
    """


    .. plot::

        from bioconvert.bam2sam import BAM2SAM
        from bioconvert import bioconvert_data
        from easydev import TempFile

        with TempFile(suffix=".sam") as fh:
            infile = bioconvert_data("test_measles.sorted.bam")
            convert = BAM2SAM(infile, fh.name)
            convert.boxplot_benchmark()


    """
    input_ext = [".bam"]
    output_ext = [".sam"]

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile:
        :param str outfile:

        command used::

            samtools view -Sbh
        """
        super(BAM2SAM, self).__init__(infile, outfile, *args, **kargs)
        self._default_method = "samtools"

    def _method_samtools(self, *args, **kwargs):
        # -S means ignored (input format is auto-detected)
        # -h means include header in SAM output
        cmd = "samtools view -Sh {} -O SAM -o {}".format(self.infile, self.outfile)
        self.execute(cmd)

    def _method_pysam(self, *args, **kwargs):
        import pysam
        pysam.sort("-o", self.outfile, self.infile)

    def _method_sambamba(self, *args, **kwargs):
        cmd = "sambamba view {} -o {}".format(self.infile, self.outfile)
        self.execute(cmd)

