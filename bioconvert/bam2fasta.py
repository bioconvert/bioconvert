"""Convert :term:`BAM` format to :term:`fasta` file"""
from bioconvert import ConvBase


class BAM2Fasta(ConvBase):
    """Bam2Fasta converter

    Wrapper of bamtools to convert bam file to fasta file.
    """
    input_ext = ['.bam']
    output_ext = ['fasta', 'fa']

    def __init__(self, infile, outfile):
        """.. rubric:: constructor
        :param str infile:
        :param str outfile:
        library used::
            pysam (samtools)
        """
        super().__init__(infile, outfile)
        self._default_method = "bamtools"

    def _method_bamtools(self, *args, **kwargs):
        """

        .. note:: fastq are split on several lines (80 characters)

        """
        # Another idea is to use pysam.bam2fq but it fails with unknown error
        #pysam.bam2fq(self.infile, save_stdout=self.outfile)
        #cmd = "samtools fastq %s >%s" % (self.infile, self.outfile)
        #self.execute(cmd)
        # !!!!!!!!!!!!!!!!!! pysam.bam2fq, samtools fastq and bamtools convert
        # give differnt answers...

        cmd = "bamtools convert -format fasta -in {0} -out {1}".format(
            self.infile, self.outfile
        )
        self.execute(cmd)

    def _method_samtools(self, *args, **kwargs):
        """
        do the conversion :term`BAM` -> :term:'Fastq` using samtools

        :return: the standard output
        :rtype: :class:`io.StringIO` object.

        .. note:: fastq are one on line
        """
        cmd = "samtools fasta {} > {}".format(self.infile, self.outfile)
        self.execute(cmd)

