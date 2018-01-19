import os
import colorlog
from Bio import SeqIO

from bioconvert import ConvBase, extensions

_log = colorlog.getLogger(__name__)


class PHYLIP2CLUSTAL(ConvBase):
    """
    Converts a sequence alignment from :term:`PHYLIP` format to :term:`CLUSTAL` format::

        converter = PHYLIP2CLUSTAL(infile, outfile)
        converter(method='biopython')

    default method = biopython
    available methods = biopython, squizz
    """

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`PHYLIP` file.
        :param str outfile: (optional) output :term:`FASTA` file
        """
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
