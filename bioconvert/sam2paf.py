# -*- coding: utf-8 -*-

###########################################################################
# Bioconvert is a project to facilitate the interconversion               #
# of life science data from one format to another.                        #
#                                                                         #
# Authors: see CONTRIBUTORS.rst                                           #
# Copyright © 2018  Institut Pasteur, Paris and CNRS.                     #
# See the COPYRIGHT file for details                                      #
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
###########################################################################

""""Convert :term:`CRAM` to :term:`BAM` format"""
import re
import os

from bioconvert import ConvBase
from easydev.multicore import cpu_count

import colorlog

from bioconvert.core.decorators import requires_nothing, requires

logger = colorlog.getLogger(__name__)


class SAM2PAF(ConvBase):
    """Convert :term:`SAM` file to :term:`PAF` file

    The :term:`SAM` and :term:`PAF` formats are described in the :ref:`formats`
    section.

    Description:

    The header of the SAM file (lines starting with @) are dropped. However, the
    length of the target is retrieved from the @SQ line that must be present.


    Consider this SAM file with two alignements only. One is aligned on the
    target (first) while the other is not (indicated by the ``*`` characters)::

        @SQ	SN:ENA|K01711|K01711.1	LN:15894
        @PG	ID:minimap2	PN:minimap2	VN:2.5-r572	CL:minimap2 -a measles.fa Hm2_GTGAAA_L005_R1_001.fastq.gz
        HISEQ:426:C5T65ACXX:5:2302:1943:2127	0	ENA|K01711|K01711.1	448	60	101M	*	00	CTTACCTTCGCATCAAGAGGTACCAACATGGAGGATGAGGCGGACCAATACTTTTCACATGATGATCCAATTAGTAGTGATCAATCCAGGTTCGGATGGTT	BCCFFFFFHHHHHIIJJJJJJIIJJJJJJJJFHIHIJJJIJIIIIGHFFFFFFEEEEEEEDDDDDFDDDDDDDDD>CDDEDEEDDDDDDCCDDDDDDDDCD	NM:i:0	ms:i:202	AS:i:202	nn:i:0	tp:A:P	cm:i:14	s1:i:94	s2:i:0
        HISEQ:426:C5T65ACXX:5:2302:4953:2090	4	*	0	0	*	*	0	0	AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAAATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAACAACCAAAAAGAGACGAACAA	CCCFFDDFAFFBHJHGGGIHIJBGGHIIJJJJJJJHGEIJGIFIIIHCBGHIJIIIIIJJHHHHEF@D@=;=,0)0&5&))+(((+((((&+(((()&&)(


    The equivalent PAF file is ::

        HISEQ:426:C5T65ACXX:5:2302:1943:2127	101	0	101	+	ENA|K01711|K01711.1	15894	447	548	101	101	60	NM:i:0	ms:i:202	AS:i:202	nn:i:0	tp:A:P	cm:i:14	s1:i:94	s2:i:0	cg:Z:101M

    In brief, the sequences are dropped. The final file is therefore smaller.
    Extra fields (starting from NM:i:0) can be dropped or kept using the
    keep_extra_field argument. Alignement with ``*`` characters are dropped.
    The first line (@SQ) is used to retrieve the length of the contigs that is
    stored in the PAF file (column 6).


    The 12 compulsary PAF fields are:

    ====== ======== ===========================================
    Col     Type    Description
    ====== ======== ===========================================
    1      string   Query sequence name
    2       int     Query sequence length
    3       int     Query start (0-based)
    4       int     Query end (0-based)
    5       char    Relative strand: "+" or "-"
    6      string   Target sequence name
    7       int     Target sequence length
    8       int     Target start on original strand (0-based)
    9       int     Target end on original strand (0-based)
    10      int     Number of residue matches
    11      int     Alignment block length
    12      int     Mapping quality (0-255; 255 for missing)
    ====== ======== ===========================================


    For developesr:

    Get the measles data from Sequana library (2 paired fastq files)::

        minimap2 measles.fa R1.fastq > approx-mapping.paf

    You can ask minimap2 to generate CIGAR at the cg tag of PAF with:

        minimap2 -c measles.fa R1.fastq > alignment.paf

    or to output alignments in the SAM format::

        minimap2 -a measles.fa R1.fastq > alignment.sam


    The SAM lines must contains 11 positional element and the NM:i and nn:i
    fields (see example above).

    """
    _default_method = "python"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input SAM file
        :param str outfile: output PAF filename

        :reference: This function is a direct translation of 
            https://github.com/lh3/miniasm/blob/master/misc/sam2paf.js (Dec. 
            2017).

        """
        super(SAM2PAF, self).__init__(infile, outfile, *args, **kargs)

    @requires_nothing
    def _method_python(self, *args, **kwargs):
        pattern = r"(\d+)([MIDSHNX=])"

        extra_fields = kwargs.get("extra_fields", "SAM")
        # TODO: what is this ?
        pri_only = kwargs.get("pri_only", True)

        skipped = 0
        reference_lengths = {}
        with open(self.infile, "r") as fin:
            with open(self.outfile, "w") as fout:
                for lineno, line in enumerate(fin.readlines()):

                    if line.startswith("@"):
                        if line.startswith("@SQ"):
                            match = re.findall(r"\tSN:(\S+)", line)
                            name = match[0] if len(match) else "unknown_reference"

                            match = re.findall(r"\tLN:(\d+)", line)
                            if len(match) == 1:
                                reference_lengths[name] = int(match[0])
                            else:
                                raise ValueError(
                                    "Could not parse SQ line to extract the length "
                                    "(LN: field missing maybe ?)")
                        continue

                    t = line.split()
                    flag = int(t[1])

                    if (t[9] != "*" and t[10] != "*" and len(t[9]) != len(t[10])):
                        raise ValueError("ERROR at line " + str(lineno) +
                            ":inconsistent SEQ and QUAL lengths - " +
                            str( len(t[9])) + " != " + str(len(t[10])))

                    if (t[2] == '*' or (flag & 4)):
                        skipped += 1
                        continue

                    # if flag is 256 and pri_only, we skip the alignment
                    if (pri_only and (flag & 0x100)):
                        continue

                    # Get the reference length for this alignment
                    if t[2] in reference_lengths:
                        tlen = reference_lengths[t[2]]
                    else:
                        raise KeyError("can't find the length of contig {}".format(t[2]))

                    # The reference is known but the length is not
                    if (tlen == -1) :
                        raise ValueError("ERROR at line " + str(lineno) + ": can't find the length of contig " + str(t[2]))

                    # TODO explain what are the nn and NM tags
                    match = re.findall(r"\tnn:i:(\d+)", line)
                    if match:
                        nn = int(match[0])
                    else:
                        nn = 0

                    match = re.findall(r"\tNM:i:(\d+)", line)
                    if match:
                        NM = int(match[0])
                        have_NM = True
                    else:
                        NM = 0
                        have_NM = False
                    NM += nn

                    # See sequana.cigar for more information
                    clip = [0, 0]
                    I = [0, 0]      # Insertion
                    D = [0, 0]      # Deletion
                    M, N = 0, 0     # Matches
                    ql, tl, mm = 0, 0, 0,
                    ext_cigar = False
                    n_cigar = 0

                    Zacc = ""
                    for count, letter in re.findall(pattern, t[5]):
                        l = int(count)
                        if (letter == 'M'):
                            M += l
                            ql += l
                            tl += l
                            ext_cigar = False
                            Zacc += "{}M".format(count)
                        elif (letter == 'I'):
                            I[0] += 1
                            I[1] += l
                            ql += l
                            Zacc += "{}I".format(count)
                        elif (letter == 'D'):
                            D[0] += 1
                            D[1] += l
                            tl += l
                            Zacc += "{}D".format(count)
                        elif (letter == 'N'):
                            N += l
                            tl += l
                        elif (letter == 'S'):
                            clip[0 if M==0 else 1] = l
                            ql += l
                        elif (letter == 'H'):
                            clip[0 if M == 0 else 1] = l
                        elif (letter == '='):
                            M += l
                            ql += l
                            tl += l
                            ext_cigar = True
                        elif (letter == 'X'):
                            M += l
                            ql += l
                            tl += l
                            mm += l
                            ext_cigar = True
                        n_cigar += 1

                    prefix_msg = "at line {}: ".format(lineno) + "{}"

                    if (n_cigar > 65535):
                         logger.warning(prefix_msg.format(str(n_cigar) +
                                        " CIGAR operations"))

                    if (tl + int(t[3]) - 1 > tlen):
                        logger.warning(prefix_msg.format("alignment end "
                            "position larger than ref length; skipped"))
                        continue

                    if (t[9] != '*' and len(t[9]) != ql) :
                        logger.warning(prefix_msg.format(
                            " SEQ length inconsistent with CIGAR(" +
                            str(len(t[9])) + " != " + str(ql) + "); skipped"))
                        continue

                    if (have_NM is False or ext_cigar):
                         NM = I[1] + D[1] + mm

                    if (NM < I[1] + D[1] + mm):
                        logger.warning(prefix_msg.format(" NM is less than the total number of gaps (" 
                            + str(NM) + " < " + str(I[1]+D[1]+mm) + ")"))
                        NM = I[1] + D[1] + mm

                    # extra information to store in the PAF after the 12
                    # extra field from original code. The insert and
                    # deletions (io/in and do/di) and mm is the number of other
                    # substitutions ? NM -I[1] -D[1]
                    extra = [
                            "mm:i:"+ str(NM-I[1]-D[1]),
                            "io:i:"+str(I[0]),
                            "in:i:"+str(I[1]),
                            "do:i:"+str(D[0]),
                            "dn:i:"+ str(D[1])]

                    match = M - (NM - I[1] - D[1])
                    blen = M + I[1] + D[1]
                    qlen = M + I[1] + clip[0] + clip[1]

                    # What does flag 16 means ?
                    if (flag & 16):
                        qs = clip[1]
                        qe = qlen - clip[0]
                    else:
                        qs = clip[0]
                        qe = qlen - clip[1]

                    ts = int(t[3]) - 1
                    te = ts + M + D[1] + N

                    ## WARNING: difference with sam2paf.js : we add and substract nn
                    ## from match and blen to agree with the output SAM file
                    ## generated by minimap2

                    # The 12 compulsary fields to have a valid PAF format
                    a = [t[0], qlen, qs, qe, "-" if flag & 16 else '+', t[2],
                         tlen, ts, te, match+nn, blen-nn, t[4]]
                    # cast to string and save in file
                    a = [str(x) for x in a]

                    # What extra fields do we want to add ?
                    # original fields found in the SAM file ?
                    if extra_fields == "SAM" and len(t)>11:
                        fout.write("\t".join(a) + "\t" + "\t".join(t[11:]) + "\tcg:Z:{}\n".format(Zacc))
                    elif extra_fields == "summary":
                        fout.write("\t".join(a) + "\t" + "\t".join(extra) + "\n")
                    elif extra_fields is None:
                        fout.write("\t".join(a) + "\n")

        self.skipped = skipped

    #@requires(external_binaries=["k8", "paftools"])
    #def _method_paftools(self, *args, **kwargs):
    #    cmd = "paftools sam2paf {} > {}".format(self.infile, self.outfile)
    #    self.execute(cmd) 
