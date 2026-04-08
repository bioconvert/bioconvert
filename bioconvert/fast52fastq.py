###########################################################################
# Bioconvert is a project to facilitate the interconversion               #
# of life science data from one format to another.                        #
#                                                                         #
# Copyright © 2018-2022  Institut Pasteur, Paris and CNRS.                #
#                                                                         #
# bioconvert is free software: you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or       #
# (at your option) any later version.                                     #
#                                                                         #
# bioconvert is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of          #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
# GNU General Public License for more details.                            #
#                                                                         #
# You should have received a copy of the GNU General Public License       #
# along with this program (COPYING file).                                 #
# If not, see <http://www.gnu.org/licenses/>.                             #
#                                                                         #
# Repository: https://github.com/bioconvert/bioconvert                    #
# Documentation: http://bioconvert.readthedocs.io                         #
###########################################################################
"""Convert :term:`FAST5` to :term:`FASTQ` format"""
import colorlog

from bioconvert import ConvBase
from bioconvert.core.decorators import requires, requires_nothing

logger = colorlog.getLogger(__name__)

__all__ = ["FAST52FASTQ"]


class FAST52FASTQ(ConvBase):
    """Convert :term:`FAST5` to :term:`FASTQ`

    FAST5 is an Oxford Nanopore Technologies (ONT) file format based on HDF5.
    When a FAST5 file contains basecalled sequence data (stored under
    ``Analyses/Basecall_1D_*/BaseCalled_template/Fastq``), this converter
    extracts those sequences and writes them as a standard FASTQ file.

    """

    _default_method = "h5py"
    _loss = False

    def __init__(self, infile, outfile):
        """
        :param str infile: The path to the input FAST5 file.
        :param str outfile: The path to the output FASTQ file.
        """
        super(FAST52FASTQ, self).__init__(infile, outfile)

    @requires(python_library="h5py")
    def _method_h5py(self, *args, **kwargs):
        """Convert :term:`FAST5` to :term:`FASTQ` using :term:`h5py`.

        Reads basecalled sequence data embedded in the FAST5 file under
        ``Analyses/Basecall_1D_*/BaseCalled_template/Fastq`` and writes
        them to a FASTQ output file.

        :raises ValueError: if no basecalled data is found in the FAST5 file.
        """
        import h5py

        with h5py.File(self.infile, "r") as f5, open(self.outfile, "w") as fq_out:
            reads_written = 0
            for read_name in f5.keys():
                read = f5[read_name]
                if "Analyses" not in read:
                    continue
                analyses = read["Analyses"]
                # Collect all Basecall_1D groups, sorted for deterministic order
                basecall_groups = sorted(
                    key for key in analyses.keys() if key.startswith("Basecall_1D")
                )
                if not basecall_groups:
                    continue
                # Use the latest basecall group
                latest = basecall_groups[-1]
                basecall = analyses[latest]
                if "BaseCalled_template" not in basecall:
                    continue
                template = basecall["BaseCalled_template"]
                if "Fastq" not in template:
                    continue
                fastq_data = template["Fastq"][()].decode("utf-8")
                # Ensure it ends with a newline
                if not fastq_data.endswith("\n"):
                    fastq_data += "\n"
                fq_out.write(fastq_data)
                reads_written += 1

        if reads_written == 0:
            raise ValueError(
                f"No basecalled FASTQ data found in {self.infile}. "
                "Please run a basecaller (e.g. Guppy) on the FAST5 file first."
            )
