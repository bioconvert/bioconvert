import os
import colorlog
from Bio import SeqIO

from bioconvert import ConvBase

_log = colorlog.getLogger(__name__)

def generate_outfile_name(infile, out_extension):
    """
    Replaces the file extension with the given one.
    :param infile: Input file
    :param out_extension: Desired extension
    :return: The filepath with the given extension
    """
    return '%s.%s' % (os.path.splitext(infile)[0], out_extension)


class FASTA2PHYLIP(ConvBase):
    """Converts a sequence alignment in :term:`FASTA` format to :term:`PHYLIP` format

    Conversion is based on Bio Python modules

    """

    input_ext = ['fa', 'fst', 'fasta', 'fn']
    output_ext = ['phylip', 'phy']

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`FASTA` file.
        :param str outfile: (optional) output :term:`PHYLIP` file
        """
        if not outfile:
            outfile = generate_outfile_name(infile, 'phylip')
        super().__init__(infile, outfile)
        self.alphabet = alphabet
        self._default_method = 'biopython'

    def _method_biopython(self):
        sequences = list(SeqIO.parse(self.infile, "fasta", alphabet=self.alphabet))
        count = SeqIO.write(sequences, self.outfile, "phylip")
        _log.debug("Converted %d records to phylip" % count)


class PHYLIP2FASTA(ConvBase):
    """Converts a sequence alignment in :term:`PHYLIP` format to :term:`FASTA` format

    Conversion is based on Bio Python modules

    """

    output_ext = ['fa', 'fst', 'fasta', 'fn']
    input_ext = ['phylip', 'phy']

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
        
    def _method_biopython(self):
        sequences = list(SeqIO.parse(self.infile, "phylip", alphabet=self.alphabet))
        count = SeqIO.write(sequences, self.outfile, "fasta")
        _log.debug("Converted %d records to fasta" % count)

