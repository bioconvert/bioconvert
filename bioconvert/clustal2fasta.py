import os
import colorlog
from Bio import SeqIO

from bioconvert import ConvBase, extensions

_log = colorlog.getLogger(__name__)


class CLUSTAL2FASTA(ConvBase):
    """
    Converts a sequence alignment from :term:`CLUSTAL` format to :term:`FASTA` format. ::

        converter = CLUSTAL2FASTA(infile, outfile)
        converter(method='biopython')

    default method = biopython
    available methods = biopython, squizz
    """

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`CLUSTAL` file.
        :param str outfile: (optional) output :term:`FASTA` file
        """
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


