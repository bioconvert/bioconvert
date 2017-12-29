"""Convert :term:`CRAM` file to :term:`BAM` file"""
import os
from bioconvert import ConvBase
from easydev.multicore import cpu_count

import colorlog
logger = colorlog.getLogger(__name__)


class SAM2PAF(ConvBase):
    """Convert :term:`SAM` file to :term:`PAF` file

    """
    input_ext = [".sam"]
    output_ext = ".paf"

    def __init__(self, infile, outfile, reference=None, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input SAM file
        :param str outfile: output PAF filename

        :reference: transalation of https://github.com/lh3/miniasm/blob/master/misc/sam2paf.js

        """
        super(SAM2PAF, self).__init__(infile, outfile, *args, **kargs)

        self._default_method = "python"

        raise NotImplementedError

    def _method_python(self, *args, **kwargs):

        pass
        """
@SQ SN:ENA|K01711|K01711.1  LN:15894
@PG ID:minimap2 PN:minimap2 VN:2.5-r572 CL:minimap2 -a measles.fa Hm2_GTGAAA_L005_R1_001.fastq.gz
HISEQ:426:C5T65ACXX:5:2302:1943:2127    0   ENA|K01711|K01711.1 448 60  101M *
0     0 CTTACCTTCGCATCAAGAGGTACCAACATGGAGGATGAGGCGGACCAATACTTTTCACATGATGATCCAATTAGTAGTGATCAATCCAGGTTCGGATGGTT BCCFFFFFHHHHHIIJJJJJJIIJJJJJJJJFHIHIJJJIJIIIIGHFFFFFFEEEEEEEDDDDDFDDDDDDDDD>CDDEDEEDDDDDDCCDDDDDDDDCD    NM:i:0  ms:i:202    AS:i:202    nn:i:0  tp:A:P  cm:i:14 s1:i:94 s2:i:0
        """
        pri_only = False
        re =  "(\d+)([MIDSHNX=])"

        with open(self.infile, "r") as fin:
            with open(self.outfile, "w") as fout:
                for line in fin.readlines():

                   if line.startswith("@"):
                        if line.startswith("@SQ"):
                            #var name = (m = /\tSN:(\S+)/.exec(line)) != null? m[1] : null;
                            #var l = (m = /\tLN:(\d+)/.exec(line)) != null? parseInt(m[1]) :null;
                            #if (name != null && l != null) len[name] = l; 
                            continue

                    #else
                    t = line.split()
                    flag = int(t[1])

                    if (t[9] != "*" and t[10] != "*" and len(t[9]) != len(t[10])):
                        raise ValueError("ERROR at line " + str(i) + ":
                            inconsistent SEQ and QUAL lengths - " str( t[9].length) + " != " +
                            str(t[10].length))

                    if (t[2] == '*' || (flag & 4)):
                        continue

                    if (pri_only && (flag&0x100)) 
                        continue

                    tlen = len(t[2])

                    if (tlen == null) :
                        raise ValueError("ERROR at line " + str(i) + ": can't find the length of contig " + str(t[2]))
                    nn = (m = /\tnn:i:(\d+)/.exec(line)) != null? parseInt(m[1]) : 0;
                    NM = (m = /\tNM:i:(\d+)/.exec(line)) != null? parseInt(m[1]) : null;
                    have_NM = NM == null? false : true
                    NM += nn
                    clip = [0, 0]
                    I = [0, 0]
                    D = [0, 0]
                    M = 0
                    N = 0
                    ql, tl, mm = 0, 0, 0,
                    ext_cigar = false

                    while ((m = re.exec(t[5])) != null) :
                        l = int(m[1])
                        if (m[2] == 'M') M += l, ql += l, tl += l, ext_cigar = false
                        else if (m[2] == 'I') ++I[0], I[1] += l, ql += l
                        else if (m[2] == 'D') ++D[0], D[1] += l, tl += l
                        else if (m[2] == 'N') N += l, tl += l
                        else if (m[2] == 'S') clip[M == 0? 0 : 1] = l, ql += l
                        else if (m[2] == 'H') clip[M == 0? 0 : 1] = l
                        else if (m[2] == '=') M += l, ql += l, tl += l, ext_cigar = true
                        else if (m[2] == 'X') M += l, ql += l, tl += l, mm += l, ext_cigar =true

                    n_cigar += 1

                    if (n_cigar > 65535)
                        _log.warning("at line " + lineno + ": " + n_cigar + " CIGAR operations")
                    if (tl + parseInt(t[3]) - 1 > tlen) {
                        _log.warning("line " + lineno + ": alignment end position larger than ref length; skipped");
                        continue
                    }
                    if (t[9] != '*' && t[9].length != ql) {
                        warn("WARNING at line " + lineno + ": SEQ length inconsistent with CIGAR
(" + t[9].length + " != " + ql + "); skipped");
                        continue;
                    }
                    if (!have_NM || ext_cigar) NM = I[1] + D[1] + mm;
                    if (NM < I[1] + D[1] + mm) {
                        warn("WARNING at line " + lineno + ": NM is less than the total number
                            of gaps (" + NM + " < " + (I[1]+D[1]+mm) + ")");
                        NM = I[1] + D[1] + mm;
                    }
                    var extra = ["mm:i:"+(NM-I[1]-D[1]), "io:i:"+I[0], "in:i:"+I[1],
"do:i:"+D[0], "dn:i:"+D[1]];
                    var match = M - (NM - I[1] - D[1]);
                    var blen = M + I[1] + D[1];
                var qlen = M + I[1] + clip[0] + clip[1];
                var qs, qe;
                if (flag&16) qs = clip[1], qe = qlen - clip[0];
                else qs = clip[0], qe = qlen - clip[1];
                var ts = parseInt(t[3]) - 1, te = ts + M + D[1] + N;
                var a = [t[0], qlen, qs, qe, flag&16? '-' : '+', t[2], tlen, ts, te, match, blen, t[4]];
                print(a.join("\t"), extra.join("\t"));
        """

