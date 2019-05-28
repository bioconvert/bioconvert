###########################################################################
# Bioconvert is a project to facilitate the interconversion               #
# of life science data from one format to another.                        #
#                                                                         #
# Authors: see CONTRIBUTORS.rst                                           #
# Copyright © 2018-2019  Institut Pasteur, Paris and CNRS.                #
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
from itertools import groupby
import math

class MAFLine(object):
    """A reader for :term:`MAF` format.

    mode refname start algsize strand refsize alignment

    ::

        a
        s ref    100 10 + 100000 ---AGC-CAT-CATT
        s contig 0   10 + 10     ---AGC-CAT-CATT

        a
        s ref    100 12 + 100000 ---AGC-CAT-CATTTT
        s contig 0   12 + 12     ---AGC-CAT-CATTTT

    The alignments are stored by pair, one item for the reference, one for the
    query. The query (second line) starts at zero.

    """
    def __init__(self, line):
        self._items = line.split()
        assert self._items[0] in 'aspq', self._items

    @property
    def name(self):
        return self._items[1]

    @property
    def sequence_size(self):
        return int(self._items[5])

    @property
    def mode(self):
        return self._items[0]

    @property
    def strand(self):
        return self._items[4]

    @property
    def alignment_start(self): # alignment start
        return int(self._items[2])

    @property
    def alignment_size(self):
        return int(self._items[3])

    @property
    def alignment(self):
        return self._items[6]

    def get_alignment_length(self):
        return len(self._items[6].replace("-", ""))


class MAF(object):
    """A reader for :term:`MAF` format."""
    def __init__(self, filename, outfile=None):
        self.filename = filename
        self.outfile = outfile

    def count_insertions(self, alnString):
        """return length without insertion, forward and reverse shift"""
        gaps = alnString.count("-")
        forwardFrameshifts = alnString.count("\\")
        reverseFrameshifts = alnString.count("/")
        letters = len(alnString) - gaps - forwardFrameshifts - reverseFrameshifts
        return letters, forwardFrameshifts, reverseFrameshifts

#    @SQ SN:NC_002929    LN:4086189
#    @PG ID:bioconvert VN:?? CL:bioconvert input.maf output.sam

    def to_sam(self):
        # identifier flag ref start qual cigar * 0 0 sequence_ref qual NM: MD:
        # AS XS RG ....

        fout = open(self.outfile, "w")
        fout.write("@HD\tVN:1.3\tSO:unsorted\n")

        msg = "maf2paf found q starting a new line."
        with open(self.filename, "r") as fin:
            # scan once to get sequence name and length.
            # only reference is of interest
            names = set()
            top = True
            for line in fin:
                if len(line.strip()) == 0:
                    top = True
                    continue

                m = MAFLine(line)
                if m.mode == "s" and top is True:
                    if m.name not in names:
                        names.add(m.name)
                        fout.write("@SQ\tSN:{}\tLN:{}\n".format(m.name,
                                                                m.sequence_size))
                    top = False

        from bioconvert import version
        fout.write("@PG\tID:{0}\tPN:{0}\tVN:{1}\tCL:{0} {2} {3}\n".format(
            "bioconvert", version, self.filename, self.outfile))

        with open(self.filename, "r") as fin:
            for line in fin:
                # skipping empty lines
                if line.strip() == "":
                    continue

                if line[0] == "a":
                    s = []
                    if "=" in line:
                        tags = dict(i.split("=") for i in line[1].split())
                    else:
                        tags = {}
                elif line[0] == "s":
                    s.append(line)
                elif line[0] == "q":
                    # quality ?
                    raise NotImplementedError(msg)
                elif line[0] == "p":
                    raise NotImplementedError(msg)
                elif this[0] == "#":
                    print(this)

                # now that we have the two lines, save into SAM file
                if len(s)>2:
                    raise NotImplementedError("mutliple alignment not implemented yet")

                if len(s) == 2:
                    ref = MAFLine(s[0])
                    query = MAFLine(s[1])

                    if ref.strand != "+":
                        raise Exception("for SAM, the 1st strand in each alignment must be +")

                    flag = self.get_flag(ref.alignment, query.strand)

                    # This is the slowest part of the code
                    cigar = get_cigar(ref, query)

                    qual = "*"

                    if "mismap" in tags:
                        mapq = tags["mismap"]
                    else:
                        mapq = "255"  # missing  254 is maximum

                    pos = int(ref.alignment_start) + 1 # convert to 1-based coordinate
                    data = [query.name, flag, ref.name, pos, mapq, cigar, "*",0,0,
                            query.alignment.replace("-","").upper(), qual]
                    #MD
                    #XS
                    #RG:Z:1
                    if "score" in tags:
                        data.append("AS:i:{}".format(tags['score']))

                    if "expect" in tags:
                        data.append("EZ:Z:{}".format(tags['expect']))


                    #try: mapq = mapqFromProb(maf.namesAndValues["mismap"])
                    #except KeyError: mapq = mapqMissing
                    editDistance = sum(1 for x, y in zip(ref.alignment, query.alignment) 
                        if x != y)
                    # no special treatment of ambiguous bases: might be a minor bug
                    editDistance = "NM:i:" + str(editDistance)

                    data.append(editDistance)

                    fout.write("\t".join([str(x) for x in data]) + "\n")
                    s = []

        if len(s):
            raise ValueError("Your MAf file seems truncated. ")

        fout.close()

    def get_flag(self, qName, query_strand):
        if qName.endswith("/1"):
            qName = qName[:-2]
            if query_strand == "+": flag = "99"  # 1 + 2 + 32 + 64
            else:              flag = "83"  # 1 + 2 + 16 + 64
        elif qName.endswith("/2"):
            qName = qName[:-2]
            if query_strand == "+": flag = "163"  # 1 + 2 + 32 + 128
            else:              flag = "147"  # 1 + 2 + 16 + 128
        else:
            if query_strand == "+": flag = "0"
            else:              flag = "16"
        return flag


def mapqFromProb(probString):
    mapqMaximum = 100
    try: p = float(probString)
    except ValueError: raise Exception("bad probability: " + probString)
    if p < 0 or p > 1: raise Exception("bad probability: " + probString)
    if p == 0: return mapqMaximum
    phred = -10 * math.log(p, 10)
    if phred >= mapqMaximum: return str(mapqMaximum)
    return str(int(round(phred)))


def cigarCategory(alignmentColumn):
    x, y = alignmentColumn
    if x == "-":
        if y == "-":
            return "P"
        else:
            return "I"
    else:
        if y == "-":
            return "D"
        elif x != y:
            return "X"
        else:
            return "M"


def cigarParts(beg, alignmentColumns, end):
    if beg:
        yield str(beg) + "H"
    # (doesn't handle translated alignments)
    for k, v in groupby(alignmentColumns, cigarCategory):
        yield str(sum(1 for _ in v)) + k
    if end:
        yield str(end) + "H"


def get_cigar(m1, m2):
    qRevStart = m2.sequence_size - m2.alignment_start - m2.alignment_size
    cigar = "".join(cigarParts(m2.alignment_start,
                               zip(m1.alignment, m2.alignment),
                               qRevStart))
    return cigar





