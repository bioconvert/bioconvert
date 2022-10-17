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
""" Convert :term:`FASTA` (scaffold) to :term:`FASTA` (contig) and :term:`AGP` formats"""
from math import log10

from bioconvert import ConvBase, requires
from bioconvert.core.base import ConvArg
from bioconvert.core.decorators import requires, requires_nothing

__all__ = ["FASTA2FASTA_AGP"]


class FASTA2FASTA_AGP(ConvBase):
    """Convert :term:`FASTA` file of scaffolds to a FASTA file of
        contigs and an :term:`AGP` file

        Method implemented in Python by bioconvert developers.

        Possible option for the future version:

        - scaffolds shorter than this length will be excluded (200)
        - scaftigs shorter than this length will be masked with "N"s (50)
        - minimum length of scaffolding stretch of Ns (default 10) to split
          scaffold into several contigs

        if input sequence is on several lines, the output contig file
        will save the sequence on a single line

        converts to upper cases

        version 2.0 (columns 9 not empty)
        columns are:


        * Column 1 contains the object, that is the identifier for the object
          assembled. This can be a chromosome, scaffold or contig. If an accession
          version number is not used, convention is to preceded a chromosome with *chr*
          (e.g. chr1) and linkage group numbers with *LG*. Contigs and scaffolds may have
          identifier that is unique within the assembly.
        * Columns 2 contains the *object_beg*: The starting coordinates of the
          component/gap on the object in column 1. These are the location in the object’s
          coordinate system, not the component’s.
        * Columns 3 is the *object_end*: The ending coordinates of the component/gap
          on the object in column 1. These are the location in the object’s coordinate
          system, not the component’s.
        * Columns 4 is the *part_number*: The line count for the components/gaps
          that make up the object described in column 1.
        * Columns 5 component_type: The sequencing status of the component. These
          typically correspond to keywords in the International Sequence Database
          (GenBank/EMBL/DDBJ) submission. Some acceptable values used here: W (WGS
            contig), N (gap with specified size), U (gap of unknown size)
        * columns 6a: component_id
        * columns 6b: gap_length
        * columns 7a: component_beg:
        * columns 7b: gap_type: if columns 5 in {N,U} specifices gap type:
            - scaffold
            - contig
            - centromere
            - short_arm
            - telomere
            - repeat
            - contamination
        * columns 8a: component_end
        * columns 8b: linkage (if column equal to N or U): tells if there is
          evidence of linkage between the adjacent lines (yes or no)
        * Columns 9a orientation: if column 5 not in {N, U}, specifies orientation
          of object in column1: + (plus), - (minus), ? (unknown), 0 (zero), na
          (irrelevant). ?, 0 and na are deprecated. use + instead (default)
        * Columns 9b linkage evidence: if column 5 is in {N,U}, linkage accepted
          values can be diverse:
            - paired-ends (paired sequences from the two ends of
              a DNA fragment, mate-pairs and molecular-barcoding.),
            -

    https://github.com/bcgsc/abyss/blob/master/bin/abyss-fatoagp

    """

    #: Default value
    _default_method = "python"
    min_scaffold_length = 200
    min_stretch_of_Ns = 10

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input FASTA file
        :param str outfile: output AGP filename

        """
        super(FASTA2FASTA_AGP, self).__init__(infile, outfile)
        self.outfile_fasta = outfile[0]
        self.outfile_agp = outfile[1]

    def _mask_scatigs(self, x, min_scatigs=10):
        # scaftigs shorter than this length will be masked with "N"s
        # AAAAAAAAAAAAA-NNNNN-GG-NNNNNN-AAAAAAAAAAAAA
        # Here the GG will be masked.
        def tt():
            if len(x) == 0:
                return "N"
            elif len(x) <= min_scatigs:
                return len(x) * "N"
            else:
                return x

        return "".join([tt(x) for x in f.split("N")])

    @requires_nothing
    def _method_python(self, *args, **kwargs):
        """Converts the input FASTA (scaffold) into FASTA (contigs) and AGP. Internal method"""
        min_scaffold_length = kwargs.get("min_scaffold_length", self.min_scaffold_length)
        stretch_of_Ns = kwargs.get("min_stretch_of_Ns", self.min_stretch_of_Ns)
        stretch_of_Ns = "N" * stretch_of_Ns

        # First, we need to figure out number of sequences
        counter = 0
        with open(self.infile, "r") as fin:
            for line in fin.readlines():
                if line.startswith(">"):
                    counter += 1

        # a simple lambda function
        ZFILL = int(log10(counter)) + 1

        def _frmt(counter):
            return "{}".format(str(counter).zfill(ZFILL))

        # We scan the input scaffold file and create (1) the contig fasta file
        # and (2) the AGP file.
        scaffold_counter = 0
        contig_counter = 0

        # given a scaffold AAANNNCCCNNNTTT :
        # This is made of 3 contigs (AAA, CCC and TTT)
        # The 2 stretch of Ns are scaffold. The linkage is paired-ends

        # NNNN should be trimmed

        # mask scaftigs shorter than -S threshold with "N"s

        with open(self.infile, "r") as fin:
            fout_fasta = open(self.outfile_fasta, "w")
            fout_agp = open(self.outfile_agp, "w")
            fout_agp.write("##agp-version   2.0\n")

            for line in fin.readlines():
                if line.startswith(">") and scaffold_counter != 0:
                    # a new line and we already read the first sequence.
                    # So, so we need to save the contig now and the
                    # AGP information into a file
                    fout_fasta.write(">contig_" + _frmt(scaffold_counter) + "\n")
                    fout_fasta.write("{}\n".format(current))

                    # Save the previous AGP line if required
                    L = len(current)
                    data = [
                        "scaffold_" + _frmt(scaffold_counter),
                        1,
                        L,
                        1,
                        "W",
                        "contig_" + _frmt(contig_counter),
                        1,
                        L,
                        "+",
                    ]
                    data = [str(x) for x in data]
                    fout_agp.write("\t".join(data) + "\n")

                    # reset the current sequence
                    current = ""

                    # now we increment the counter
                    scaffold_counter += 1

                elif line.startswith(">") is False:
                    # accumulating the sequence with upper case
                    current += line.strip().upper()
                    # Do we have Ns
                    indexN = current.find(stretch_of_Ns)
                    if indexN != -1:
                        print("Found Ns")

                elif line.startswith(">"):
                    # This is the first line of the file
                    current = ""
                    scaffold_counter += 1
                else:
                    raise IOError("Not a valid FASTA file")

            # finally, the last line in the AGP file
            fout_fasta.write(">contig_" + _frmt(scaffold_counter) + "\n")
            fout_fasta.write("{}\n".format(current))
            L = len(current)
            data = [
                "scaffold_" + _frmt(scaffold_counter),
                1,
                L,
                1,
                "W",
                "contig_" + _frmt(contig_counter),
                1,
                L,
                "+",
            ]
            data = [str(x) for x in data]
            fout_agp.write("\t".join(data) + "\n")
            fout_agp.close()
            fout_fasta.close()

    @classmethod
    def get_additional_arguments(cls):
        yield ConvArg(
            names="--min-scaffold-length",
            default=FASTA2FASTA_AGP.min_scaffold_length,
            type=int,
            help="minimum scaffold length",
        )
        yield ConvArg(
            names="--min-stretch-of-N",
            default=FASTA2FASTA_AGP.min_stretch_of_Ns,
            type=int,
            help="minimum stretch of Ns to split a scaffold",
        )

    @staticmethod
    def get_IO_arguments():
        yield ConvArg(
            names="input_file",
            default=None,
            type=ConvArg.file,
            help="Path to the input scaffold FASTA file.",
        )
        yield ConvArg(
            names="output_file",
            nargs=2,
            default=None,
            type=ConvArg.file,
            output_argument=True,
            help="contig FASTA file followed by the AGP file.",
        )
