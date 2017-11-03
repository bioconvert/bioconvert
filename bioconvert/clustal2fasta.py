import os
import colorlog
from Bio import SeqIO

from bioconvert import ConvBase, generate_outfile_name

_log = colorlog.getLogger(__name__)


class CLUSTAL2FASTA(ConvBase):
    """
    Converts a sequence alignment from :term:`CLUSTAL` format to :term:`FASTA` format. ::

        converter = CLUSTAL2FASTA(infile, outfile)
        converter(method='biopython')

    default method = biopython
    available methods = biopython, squizz
    """

    input_ext = ['clustal', 'clw']
    output_ext = ['fasta', 'fa', 'fst']

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`CLUSTAL` file.
        :param str outfile: (optional) output :term:`FASTA` file
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
        count = SeqIO.write(sequences, self.outfile, "fasta")
        _log.info("Converted %d records to fasta" % count)


    def _method_squizz(self, threads=None):
        """
        Convert :term:`CLUSTAL` file in :term:`FASTA` format using squizz tool.

        :param threads: not used.
        """
        cmd = 'squizz -c FASTA {infile} > {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)


class FASTA2CLUSTAL(ConvBase):
    """
    Converts a sequence alignment from :term:`FASTA` format to :term:`CLUSTAL` format::

        converter = FASTA2CLUSTAL(infile, outfile)
        converter(method='biopython')

    default method = biopython
    available methods = biopython, squizz
    """

    input_ext = ['fasta', 'fa', 'fst']
    output_ext = ['clustal', 'clw']

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`FASTA` file.
        :param str outfile: (optional) output :term:`CLUSTAL` file
        """
        if not outfile:
            outfile = generate_outfile_name(infile, 'clustal')
        super().__init__(infile, outfile)
        self.alphabet = alphabet
        self._default_method = 'biopython'


    def _method_biopython(self, threads=None):
        """
        Convert :term:`FASTA` interleaved file in :term:`CLUSTAL` format using biopython.

        :param threads: not used
        """
        sequences = list(SeqIO.parse(self.infile, "fasta", alphabet=self.alphabet))
        count = SeqIO.write(sequences, self.outfile, "clustal")
        _log.info("Converted %d records to clustal" % count)


    def _method_squizz(self, threads=None):
        """
        Convert :term:`FASTA` file in :term:`CLUSTAL` format using squizz tool.

        :param threads: not used
        """
        cmd = 'squizz -c CLUSTAL {infile} > {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)