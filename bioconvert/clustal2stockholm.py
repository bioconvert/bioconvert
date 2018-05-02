import os
import colorlog
from Bio import SeqIO

from bioconvert import ConvBase
from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)


class CLUSTAL2STOCKHOLM(ConvBase):
    """
    Converts a sequence alignment from :term:`CLUSTAL` format to :term:`STOCKHOLM` format. ::

        converter = CLUSTAL2STOCKHOLM(infile, outfile)
        converter(method='biopython')

    default method = biopython
    available methods = biopython, squizz
    """
    _default_method = 'biopython'

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`CLUSTAL` file.
        :param str outfile: (optional) output :term:`STOCKHOLM` file
        """
        super().__init__(infile, outfile)
        self.alphabet = alphabet

    @requires(python_library="biopython")
    def _method_biopython(self, threads=None, *args, **kwargs):
        """
        Convert :term:`CLUSTAL` interleaved file in :term:`PHYLIP` format using biopython.

        :param threads: not used.
        """
        sequences = list(SeqIO.parse(self.infile, "clustal", alphabet=self.alphabet))
        count = SeqIO.write(sequences, self.outfile, "stockholm")
        _log.info("Converted %d records to stockholm" % count)

    @requires("squizz")
    def _method_squizz(self, threads=None, *args, **kwargs):
        """
        Convert :term:`CLUSTAL` file in :term:`STOCKHOLM` format using squizz tool.

        :param threads: not used.
        """
        cmd = 'squizz -c STOCKHOLM {infile} > {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)


