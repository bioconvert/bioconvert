import os
import colorlog
from Bio import SeqIO

from bioconvert import ConvBase, generate_outfile_name, check_tool, install_tool

_log = colorlog.getLogger(__name__)


class FASTA2NEXUS(ConvBase):
    """Converts a sequence alignment in :term:`FASTA` format to :term:`NEXUS` format

    Conversion is based on Bio Python modules

    """

    input_ext = ['fa', 'fst', 'fasta', 'fn']
    output_ext = ['nexus', 'nx']

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`FASTA` file.
        :param str outfile: (optional) output :term:`NEXUS` file
        """
        if not outfile:
            outfile = generate_outfile_name(infile, 'nexus')
        super().__init__(infile, outfile)
        self.alphabet = alphabet
        self._default_method = 'goalign'

    def _method_goalign(self, threads=None, *args, **kwargs):
        """
        Convert fasta file in Nexus format using goalign tool.
        https://github.com/fredericlemoine/goalign

        The fasta file must be an alignemnt file, yhis mean all the sequences must
        have the same length (with the gap) otherwise an error will be raised
        """
        if not check_tool('goalign'):
            install_tool('goalign')
        cmd = 'goalign reformat nexus -i {infile} -o {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)
