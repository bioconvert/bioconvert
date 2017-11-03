import os
import colorlog
from Bio import SeqIO

from bioconvert import ConvBase, generate_outfile_name

_log = colorlog.getLogger(__name__)


class CLUSTAL2PHYLIP(ConvBase):
    """
    Converts a sequence alignment from :term:`CLUSTAL` format to :term:`PHYLIP` format. ::

        converter = CLUSTAL2PHYLIP(infile, outfile)
        converter(method='biopython')

    default method = biopython
    available methods = biopython, squizz
    """

    input_ext = ['clustal', 'clw']
    output_ext = ['phylip', 'phy']

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`CLUSTAL` file.
        :param str outfile: (optional) output :term:`PHYLIP` file
        """
        if not outfile:
            outfile = generate_outfile_name(infile, 'phylip')
        super().__init__(infile, outfile)
        self.alphabet = alphabet
        self._default_method = 'biopython'


    def _method_biopython(self, threads=None):
        """
        Convert :term:`CLUSTAL` interleaved file in :term:`PHYLIP` format using biopython.

        :param threads: not used.
        """
        sequences = list(SeqIO.parse(self.infile, "clustal", alphabet=self.alphabet))
        count = SeqIO.write(sequences, self.outfile, "phylip")
        _log.info("Converted %d records to phylip" % count)


    def _method_squizz(self, threads=None):
        """
        Convert :term:`CLUSTAL` interleaved file in :term:`PHYLIP` format using squizz tool.

        :param threads: not used.
        """
        cmd = 'squizz -c PHYLIPI {infile} > {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)


class PHYLIP2CLUSTAL(ConvBase):
    """
    Converts a sequence alignment from :term:`PHYLIP` format to :term:`CLUSTAL` format::

        converter = PHYLIP2CLUSTAL(infile, outfile)
        converter(method='biopython')

    default method = biopython
    available methods = biopython, squizz
    """

    input_ext = ['phylip', 'phy']
    output_ext = ['clustal', 'clw']

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`PHYLIP` file.
        :param str outfile: (optional) output :term:`FASTA` file
        """
        if not outfile:
            outfile = generate_outfile_name(infile, 'fasta')
        super().__init__(infile, outfile)
        self.alphabet = alphabet
        self._default_method = 'biopython'


    def _method_biopython(self, threads=None):
        """
        Convert :term:`PHYLIP` interleaved file in :term:`CLUSTAL` format using biopython.

        :param threads: not used
        """
        sequences = list(SeqIO.parse(self.infile, "phylip", alphabet=self.alphabet))
        count = SeqIO.write(sequences, self.outfile, "clustal")
        _log.info("Converted %d records to clustal" % count)


    def _method_squizz(self, threads=None):
        """
        Convert :term:`PHYLIP` interleaved file in :term:`CLUSTAL` format using squizz tool.

        :param threads: not used
        """
        cmd = 'squizz -c CLUSTAL {infile} > {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)