import os
import colorlog
from Bio import SeqIO

from bioconvert import ConvBase, generate_outfile_name, check_tool, install_tool

_log = colorlog.getLogger(__name__)


class NEXUS2FASTA(ConvBase):
    """
    Converts a sequence alignment from :term:`NEXUS` format to :term:`FASTA` format. ::
    """

    input_ext = ['nexus', 'nx']
    output_ext = ['fasta', 'fa']

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`NEXUS` file.
        :param str outfile: (optional) output :term:`FASTA` file
        """
        if not outfile:
            outfile = generate_outfile_name(infile, 'fasta')
        super().__init__(infile, outfile)
        self.alphabet = alphabet
        self._default_method = 'goalign'


    def _method_goalign(self, threads=None, *args, **kwargs):
        """
        Convert :term:`NEXUS` interleaved file in :term:`FASTA` format using goalign tool.
        https://github.com/fredericlemoine/goalign

        :param threads: not used.
        """
        if not check_tool('goalign'):
            install_tool('goalign')
        cmd = 'goalign reformat fasta -i {infile} -o {outfile} -x'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)
