import os
import colorlog
from Bio import SeqIO

from bioconvert import ConvBase, generate_outfile_name, check_tool, install_tool

_log = colorlog.getLogger(__name__)


class PHYLIP2NEXUS(ConvBase):
    """
    Converts a sequence alignment from :term:`PHYLIP` format to :term:`NEXUS` format. ::
    """

    input_ext = ['phylip', 'phy']
    output_ext = ['nexus', 'nx']

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`PHYLIP` file.
        :param str outfile: (optional) output :term:`NEXUS` file
        """
        if not outfile:
            outfile = generate_outfile_name(infile, 'nexus')
        super().__init__(infile, outfile)
        self.alphabet = alphabet
        self._default_method = 'goalign'


    def _method_goalign(self, threads=None, *args, **kwargs):
        """
        Convert :term:`PHYLIP` interleaved file in :term:`NEXUS` format using goalign tool.
        https://github.com/fredericlemoine/goalign

        :param threads: not used.
        """
        if not check_tool('goalign'):
            install_tool('goalign')
        cmd = 'goalign reformat nexus -i {infile} -o {outfile} -p'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)
