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
"""Convert :term:`FASTQ` to :term:`FASTA` format"""
import mmap

from mappy import fastx_read

from bioconvert import ConvBase, bioconvert_script

# from bioconvert.core.base import ConvArg
from bioconvert.core.decorators import compressor, in_gz, requires, requires_nothing


class FASTQ2FASTA(ConvBase):
    """Convert :term:`FASTQ` to :term:`FASTA`

    This converter has lots of methods. Some of them have also been
    removed or commented with time. BioPython for instance is commented
    due to poo performance compared to others. Does not mean that it is
    not to be considered. Performances are decrease due to lot of sanity checks.

    Similarly, bioawk and python_external method are commented because
    redundant with other equivalent method.

    """

    # use readfq for now because pure python are fast enough
    # for production, could use seqtk which seems the fastest method
    # though. Make sure that the default handles also the compresssion
    # input_ext = extensions.extensions.fastq
    # output_ext =  extensions.fasta
    #: default value
    _default_method = "bioconvert"
    _threading = True  # not used but required for the main benchmark

    def __init__(self, infile, outfile):
        """
        :param str infile: The path to the input FASTA file.
        :param str outfile: The path to the output file.
        """
        super(FASTQ2FASTA, self).__init__(infile, outfile)

    @staticmethod
    def just_name(record):
        """
        This method takes a Biopython sequence record *record*
        and returns its name. The comment part is not included.
        """
        return record.id

    @staticmethod
    def unwrap_fasta(infile, outfile, strip_comment=False):
        """
        This method reads fasta sequences from *infile*
        and writes them unwrapped in *outfile*. Used in the test suite.

        :param str infile: The path to the input FASTA file.
        :param str outfile: The path to the output file.
        """
        from Bio import SeqIO
        from Bio.SeqIO import FastaIO

        with open(outfile, "w") as fasta_out:
            if strip_comment:
                FastaIO.FastaWriter(fasta_out, wrap=None, record2title=FASTQ2FASTA.just_name).write_file(
                    SeqIO.parse(infile, "fasta")
                )
            else:
                FastaIO.FastaWriter(fasta_out, wrap=None).write_file(SeqIO.parse(infile, "fasta"))

    # Adapted from the readfq code by Heng Li
    # (https://raw.githubusercontent.com/lh3/readfq/master/readfq.py)
    @staticmethod
    def readfq(fp):  # this is a generator function
        last = None  # this is a buffer keeping the last unprocessed line
        while True:  # mimic closure; is it a bad idea?
            if not last:  # the first record or a record following a fastq
                for l in fp:  # search for the start of the next record
                    if l[0] == "@":  # fastq header line
                        last = l[:-1]  # save this line
                        break
            if not last:
                break
            header, seqs, last = last[1:], [], None
            for l in fp:  # read the sequence
                if l[0] in "@+":
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            seq, leng, seqs = "".join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield header, seq, "".join(seqs)  # yield a fastq record
                    break

    # @requires(python_library="biopython")
    # @compressor
    # def _method_biopython(self, *args, **kwargs):
    #     """For this method we use the biopython package Bio.SeqIO.

    #     `Bio.SeqIO Documentation <https://biopython.org/docs/1.76/api/Bio.SeqIO.html>`_"""
    #     from Bio import SeqIO
    #     records = SeqIO.parse(self.infile, 'fastq')
    #     SeqIO.write(records, self.outfile, 'fasta')

    @requires(external_binary="seqtk")
    @compressor
    def _method_seqtk(self, *args, **kwargs):
        # support gz files natively
        """We use the Seqtk library.

        `Documentation of the Seqtk method <https://github.com/lh3/seqtk>`_"""
        cmd = "seqtk seq -A {} > {}".format(self.infile, self.outfile)
        self.execute(cmd)

    @requires(external_binary="seqkit")
    @compressor
    def _method_seqkit(self, *args, **kwargs):
        # support gz files natively
        """We use the Seqkit library.

        `Documentation of the Seqkit method <https://github.com/shenwei356/seqkit>`_"""
        cmd = "seqkit fq2fa {} -j {} > {}".format(self.infile, self.threads, self.outfile)
        self.execute(cmd)

    @requires_nothing
    @compressor
    def _method_readfq(self, *args, **kwargs):
        """This method is inspired by Readfq coded by Heng Li.

        `original Readfq method <https://github.com/lh3/readfq>`_"""
        with open(self.outfile, "w") as fasta, open(self.infile, "r") as fastq:
            for (name, seq, _) in FASTQ2FASTA.readfq(fastq):
                fasta.write(">{}\n{}\n".format(name, seq))

    # Does not give access to the comment part of the header
    @requires(python_library="mappy")
    @compressor
    def _method_mappy(self, *args, **kwargs):
        """This method provides a fast and accurate C program to align genomic
        sequences and transcribe nucleotides.

        `mappy method <https://pypi.org/project/mappy/>`_"""
        with open(self.outfile, "w") as fasta:
            for (name, seq, _) in fastx_read(self.infile):
                fasta.write(">{}\n{}\n".format(name, seq))

    @requires("awk")
    @compressor
    def _method_awk(self, *args, **kwargs):
        """Here we are using the awk method.

        .. note::  Another method with awk has been tested but is less efficient. Here is which one was used:

            ::

                box.awkcmd = \"\"\"awk \'{{if(NR%4==1) {{printf(\">%s\\n\",substr($0,2));}} else if(NR%4==2) print;}}\' \"\"\"

        `awk documentation <https://www.gnu.org/software/gawk/manual/gawk.html>`_"""
        awkcmd = """awk '{{print ">"substr($0,2);getline;print;getline;getline}}'"""
        cmd = "{} {} > {}".format(awkcmd, self.infile, self.outfile)
        self.execute(cmd)

    @requires("mawk")
    @compressor
    def _method_mawk(self, *args, **kwargs):
        """This variant of the awk method uses mawk, a lighter and faster
        implementation of awk.

        .. note :: Other methods with mawk have been tested but are less efficient. Here are which ones were used:

            ::

                mawkcmd_v2 = \"\"\"mawk \'{{if(NR%4==1) {{printf(\">%s\\n\",substr($0,2));}} else if(NR%4==2) print;}}\' \"\"\"
                mawkcmd_v3 = \"\"\"mawk \'(++n<=0){next}(n!=1){print;n=-2;next}{print\">\"substr($0,2)}\'\"\"\"

        `mawk documentation <https://invisible-island.net/mawk/manpage/mawk.html>`_"""
        awkcmd = """mawk '{{print ">"substr($0,2);getline;print;getline;getline}}'"""
        cmd = "{} {} > {}".format(awkcmd, self.infile, self.outfile)
        self.execute(cmd)

    @requires("sed")
    @compressor
    def _method_sed(self, *args, **kwargs):
        """This method uses the UNIX function sed which is a non-interactive editor.

        .. note::  Another method with sed has been tested but is less efficient. Here is which one was used:

            ::

                cmd = \"\"\"sed -n \'s/^@/>/p;n;p;n;n\'\"\"\"

        `sed documentation <https://www.gnu.org/software/sed/manual/sed.html>`_"""
        cmd = """sed -n '1~4s/^@/>/p;2~4p' """
        cmd = "{} {} > {}".format(cmd, self.infile, self.outfile)
        self.execute(cmd)

    @requires("perl")
    @compressor
    def _method_perl(self, *args, **kwargs):
        """This method uses the perl command which will call the
        \"fastq2fasta.pl\" script.

        `Perl documentation <https://perldoc.perl.org/>`_"""
        perlcmd = "perl {}".format(bioconvert_script("fastq2fasta.pl"))
        cmd = "{} {} {}".format(perlcmd, self.infile, self.outfile)
        self.execute(cmd)

    @requires_nothing
    @compressor
    def _method_bioconvert(self, *args, **kwargs):
        """Bioconvert implementation in pure Python.
        This is the default method because it is the fastest."""
        with open(self.infile, "r+") as inp:

            with open(self.outfile, "wb") as out:
                mapp = mmap.mmap(inp.fileno(), 0)
                line = mapp.readline()
                while line:
                    out.write(b">")
                    out.write(line[1:])
                    out.write(mapp.readline())
                    mapp.readline()
                    mapp.readline()
                    line = mapp.readline()
                mapp.close()

    """@requires_nothing
    def _method_python_external(self, *args, **kwargs):
        pycmd = "python {}".format(bioconvert_script("fastq2fasta.py"))
        cmd = "{} {} {}".format(pycmd, self.infile, self.outfile)
        self.execute(cmd)
    """

    """@classmethod
    def get_additional_arguments(cls):
        yield ConvArg(
            names="--quality-file",
            nargs=1,
            default=None,
            type=ConvArg.file,
            output_argument=True,
            help="The path to the quality file.",
        )
    """
