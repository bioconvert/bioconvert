import os
import colorlog
from Bio import SeqIO

from bioconvert import ConvBase, generate_outfile_name, check_tool, install_tool

_log = colorlog.getLogger(__name__)


class PHYLOXML2NEXUS(ConvBase):
    """
    Converts a tree file from :term:`PHYLOXML` format to :term:`NEXUS` format. ::
    """

    input_ext = ['phyloxml', 'xml']
    output_ext = ['nexus', 'nx']

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`PHYLOXML` file.
        :param str outfile: (optional) output :term:`NEXUS` file
        """
        if not outfile:
            outfile = generate_outfile_name(infile, 'nexus')
        super().__init__(infile, outfile)
        self.alphabet = alphabet
        self._default_method = 'gotree'

    def _method_gotree(self, threads=None, *args, **kwargs):
        """
        Convert :term:`PHYLOXML`  file in :term:`NEXUS` format using gotree tool.
        https://github.com/fredericlemoine/gotree

        :param threads: not used.
        """
        if not check_tool('gotree'):
            install_tool('gotree')
        cmd = 'gotree reformat nexus -i {infile} -o {outfile} -f phyloxml'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)
