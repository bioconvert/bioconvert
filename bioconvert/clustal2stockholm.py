import os
import colorlog
from Bio import SeqIO

from bioconvert import ConvBase, generate_outfile_name

_log = colorlog.getLogger(__name__)


class CLUSTAL2STOCKHOLM(ConvBase):
    """
    Converts a sequence alignment from :term:`CLUSTAL` format to :term:`STOCKHOLM` format. ::

        converter = CLUSTAL2STOCKHOLM(infile, outfile)
        converter(method='biopython')

    default method = biopython
    available methods = biopython, squizz
    """

    input_ext = ['clustal', 'clw']
    output_ext = ['stockholm']

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`CLUSTAL` file.
        :param str outfile: (optional) output :term:`STOCKHOLM` file
        """
        if not outfile:
            outfile = generate_outfile_name(infile, 'clustal')
        super().__init__(infile, outfile)
        self.alphabet = alphabet
        self._default_method = 'biopython'


    def _method_biopython(self, threads=None):
        """
        Convert :term:`CLUSTAL` interleaved file in :term:`PHYLIP` format using biopython.

        :param threads: not used.
        """
        sequences = list(SeqIO.parse(self.infile, "clustal", alphabet=self.alphabet))
        count = SeqIO.write(sequences, self.outfile, "stockholm")
        _log.info("Converted %d records to stockholm" % count)


    def _method_squizz(self, threads=None):
        """
        Convert :term:`CLUSTAL` file in :term:`STOCKHOLM` format using squizz tool.

        :param threads: not used.
        """
        cmd = 'squizz -c STOCKHOLM {infile} > {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)


class STOCKHOLM2CLUSTAL(ConvBase):
    """
    Converts a sequence alignment from :term:`STOCKHOLM` format to :term:`CLUSTAL` format::

        converter = STOCKHOLM2CLUSTAL(infile, outfile)
        converter(method='biopython')

    default method = biopython
    available methods = biopython, squizz
    """

    input_ext = ['stockholm']
    output_ext = ['clustal', 'clw']

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`STOCKHOLM` file.
        :param str outfile: (optional) output :term:`CLUSTAL` file
        """
        if not outfile:
            outfile = generate_outfile_name(infile, 'clustal')
        super().__init__(infile, outfile)
        self.alphabet = alphabet
        self._default_method = 'biopython'


    def _method_biopython(self, threads=None):
        """
        Convert :term:`STOCKHOLM` interleaved file in :term:`CLUSTAL` format using biopython.

        :param threads: not used
        """
        sequences = list(SeqIO.parse(self.infile, "stockholm", alphabet=self.alphabet))
        count = SeqIO.write(sequences, self.outfile, "clustal")
        _log.info("Converted %d records to clustal" % count)


    def _method_squizz(self, threads=None):
        """
        Convert :term:`STOCKHOLM` file in :term:`CLUSTAL` format using squizz tool.

        :param threads: not used
        """
        cmd = 'squizz -c CLUSTAL {infile} > {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)