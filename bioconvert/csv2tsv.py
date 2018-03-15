"""Convert :term:`CSV` format to :term:`TSV` file"""

import colorlog

try:
    import pandas as pd
except:
    pass

from bioconvert import tsv2csv

logger = colorlog.getLogger(__name__)


class CSV2TSV(tsv2csv.TSV2CSV):
    """Bam2Json converter

    Convert bam file to json file.
    """
    _default_method = "python"

    def __init__(self, infile, outfile):
        """.. rubric:: constructor
        :param str infile:
        :param str outfile:
        """
        super().__init__(infile, outfile)

    def _method_python(self, in_sep=',', out_sep='\t', *args, **kwargs):
        """
        do the conversion :term`CSV` -> :term:'TSV` using standard Python modules
        """
        super()._method_python(in_sep=in_sep, out_sep=out_sep, *args, **kwargs)

    def _method_panda(self, in_sep=',', out_sep='\t', *args, **kwargs):
        """
        do the conversion :term`CSV` -> :term:'TSV` using Panda modules
        """
        super()._method_panda(in_sep=in_sep, out_sep=out_sep, *args, **kwargs)
