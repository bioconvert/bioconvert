"""Convert :term:`TSV` format to :term:`CSV` file"""
import csv

import colorlog

try:
    import pandas as pd
except:
    pass

from bioconvert import ConvBase

logger = colorlog.getLogger(__name__)


class TSV2CSV(ConvBase):
    """TSV2CSV converter

    Convert tsv file to csv file.
    """
    _default_method = "python"

    def __init__(self, infile, outfile):
        """.. rubric:: constructor
        :param str infile:
        :param str outfile:
        """
        super().__init__(infile, outfile)

    def _method_python(self, in_sep='\t', out_sep=',', line_terminator='\n', *args, **kwargs):
        """
        do the conversion :term`TSV` -> :term:'CSV` using csv module to read and write

        """
        with open(self.infile, "r") as in_stream, open(self.outfile, "w") as out_stream:
            writer = csv.writer(out_stream, delimiter=out_sep, lineterminator=line_terminator)
            reader = csv.reader(in_stream, delimiter=in_sep)
            for row in reader:
                writer.writerow(row)

    def _method_python_v2(self, in_sep='\t', out_sep=',', line_terminator='\n', *args, **kwargs):
        """
        do the conversion :term`TSV` -> :term:'CSV` using csv module to read, and writing directly into the file. Note
        that his method can't escape nor quote output char

        """
        with open(self.infile, "r") as in_stream, open(self.outfile, "w") as out_stream:
            reader = csv.reader(in_stream, delimiter=in_sep)
            for row in reader:
                out_stream.write(out_sep.join(row))
                out_stream.write(line_terminator)

    def _method_panda(self, in_sep='\t', out_sep=',', line_terminator='\n', *args, **kwargs):
        """
        do the conversion :term`TSV` -> :term:'CSV` using Panda modules

        """
        pd.read_csv(
            self.infile,
            sep=in_sep,
        ) \
            .to_csv(
            self.outfile,
            sep=out_sep,
            line_terminator=line_terminator,
            index=False,
            header='infer'
        )
