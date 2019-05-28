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
"""Sniffer for all formats included in Bioconvert"""
from bioconvert.core.extensions import extensions

import colorlog
_log = colorlog.getLogger(__name__)


class Sniffer(object):
    """Sniffer for formats included in Bioconvert

    :attr:`bioconvert.core.extensions.extensions`

    ::

        >>> from bioconvert import Sniffer
        >>> s =  Sniffer()
        >>> s.sniff("test.clustal")
        "clustal"

    TSV and CSV matchess many formats. For example a BED file is compatible with
    the TSV format.

    If a file matches several 2 formats and includes a TSV/CSV, we ignore the
    TSV/CSV.


    """
    formats = sorted(extensions.keys())

    def __init__(self):
        pass

    def sniff(self, filename):
        """Return first frmt found to be compatible with the input file"""

        candidates = []
        for frmt in self.formats:
            _log.debug("Trying {}".format(frmt))
            func = getattr(self, "is_{}".format(frmt))
            try:
                ret = func(filename)                
                if ret is True:
                    candidates.append(frmt)
            except NotImplementedError:
                pass
            #except Exception as err:
            #    raise(err)
        if "tsv" in candidates or "csv" in candidates:
            _log.warning("Ignore TSV/CSV: {}".format(candidates))
            candidates = [x for x in candidates if x not in ["tsv", "csv"]]


        if len(candidates) == 0:
            return None
        elif len(candidates) == 1:
            return candidates[0]
        else:
            _log.warning("Sniffer found several candidates: {}".format(candidates))
            return candidates

    def _is_blank_line(self, line):
        line = line.strip()
        if len(line) == 0:
            return True
        else:
            return False

    def is_abi(self, filename):
        raise NotImplementedError

    def is_bam(self, filename):
        raise NotImplementedError

    def is_bcf(self, filename):
        raise NotImplementedError

    def is_bed(self, filename):
        raise NotImplementedError

    def is_bedgraph(self, filename):
        raise NotImplementedError

    def is_bigwig(self, filename):
        raise NotImplementedError

    def is_bigbed(self, filename):
        raise NotImplementedError

    def is_bplink(self, filename):
        raise NotImplementedError

    def is_bz2(self, filename):
        raise NotImplementedError

    def is_cdao(self, filename):
        raise NotImplementedError

    def is_csv(self, filename):

        try:
            import pandas as pd
            df = pd.read_csv(filename, sep=",")
            if len(df.columns) > 1:
                return True
        except:
            pass

    def is_cram(self, filename):
        raise NotImplementedError

    def is_clustal(self, filename):
        with open(filename, "r") as fin:
            try:
                line = fin.readline().strip()
                if self._is_blank_line(line):
                    pass
                elif line.startswith("CLUSTAL"):
                    return True
            except:
                return False

    def is_dsrc(self, filename):
        raise NotImplementedError

    def is_embl(self, filename):
        raise NotImplementedError

    def is_fasta(self, filename):
        raise NotImplementedError

    def is_fastq(self, filename):
        raise NotImplementedError

    def is_genbank(self, filename):

        with open(filename, "r") as fin:
            try:
                line = fh.readline().strip()
                if self._is_blank_line(line):
                    pass
                else:
                    data = line.split()
                    data[0] in ['LOCUS']
                return True
            except:
                return False

    def is_gfa(self, filename):
        raise NotImplementedError

    def is_gff2(self, filename):
        raise NotImplementedError

    def is_gff3(self, filename):
        raise NotImplementedError

    def is_gz(self, filename):
        raise NotImplementedError

    def is_json(self, filename):
        raise NotImplementedError

    def is_maf(self, filename):
        with open(filename, "r") as fin:
            try:
                # read at most 50 lines and figure out whether
                # some lines starts with a or s
                # we get rid of the comments.
                # Read 5000 characters at most.
                data = fin.readlines(5000)
                comments = [line for line in data if line.startswith('#')]
                data = [line.strip() for line in data if line.startswith('#') is False]

                # get rid of blank lines
                data = [line for line in data if len(line.strip())!=0]
                starts = [line[0:2] for line in data]

                if len(starts) == 0: 
                    return False

                # line must start with one of i, e, q, a, s letter
                for x in starts:
                    assert x in ['a ', 's ', 'e ', 'q ', 'i ']
                return True
            except Exception as err:
                _log.debug(err)
                return False

    def is_newick(self, filename):
        raise NotImplementedError

    def is_nexus(self, filename):
        raise NotImplementedError

    def is_ods(self, filename):
        raise NotImplementedError

    def is_paf(self, filename):
        raise NotImplementedError

    def is_phylip(self, filename):

        with open(filename, "r") as fin:
            # First, we figure out the dimensions of the alignemnt.
            # we should find 2 integers
            # blank lines are forbidden in  principle between header an
            # alignment
            try:
                header = fin.readline().strip()
                first = fin.readline().strip()
                m, n = header.split()
                m = int(m)
                n = int(n)
                name, seq = first.split(" ", 1)
                seq = seq.replace(" ", "")
                # we identify each alignement and check that the length are
                # identical and equal to n
                for this in range(1, m-1): # -1 since we already read 1 line
                    nextline = fin.readline().strip()
                    name, seq2 = nextline.split(" ", 1)
                    seq2= seq2.replace(" ","")
                    assert len(seq) == len(seq2), "not same length"
                return True
            except Exception as err:
                #print(err)
                return False

    def is_phyloxml(self, filename):
        raise NotImplementedError

    def is_plink(self, filename):
        raise NotImplementedError

    def is_qual(self, filename):
        raise NotImplementedError

    def is_sam(self, filename):
        raise NotImplementedError

    def is_scf(self, filename):
        raise NotImplementedError

    def is_sra(self, filename):
        # not need. This is not a format.
        raise NotImplementedError

    def is_stockholm(self, filename):
        with open(filename, "r") as fin:
            try:
                header = fin.readline().strip()
                assert "STOCKHOLM" in header

                return True
            except:
                return False

    def is_twobit(self, filename):
        raise NotImplementedError

    def is_tsv(self, filename):
        try:
            import pandas as pd
            df = pd.read_csv(filename, sep="\s+")
            if len(df.columns) > 1:
                return True
        except:
            pass

    def is_vcf(self, filename):
        raise NotImplementedError

    def is_wiggle(self, filename):
        raise NotImplementedError

    def is_wig(self, filename):
        raise NotImplementedError

    def is_xls(self, filename):
        raise NotImplementedError

    def is_xlsx(self, filename):
        raise NotImplementedError

    def is_xmfa(self, filename):
        raise NotImplementedError

    def is_yaml(self, filename):
        raise NotImplementedError


