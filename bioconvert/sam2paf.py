"""Convert :term:`CRAM` file to :term:`BAM` file"""
import re
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

        raise NotImplementedError("Need to decide what to do with unmapped reads ! ")

    def _method_python(self, *args, **kwargs):

        """
@SQ SN:ENA|K01711|K01711.1  LN:15894
@PG ID:minimap2 PN:minimap2 VN:2.5-r572 CL:minimap2 -a measles.fa Hm2_GTGAAA_L005_R1_001.fastq.gz

HISEQ:426:C5T65ACXX:5:2302:1943:2127    0   ENA|K01711|K01711.1 448 60  101M *
0     0 CTTACCTTCGCATCAAGAGGTACCAACATGGAGGATGAGGCGGACCAATACTTTTCACATGATGATCCAATTAGTAGTGATCAATCCAGGTTCGGATGGTT BCCFFFFFHHHHHIIJJJJJJIIJJJJJJJJFHIHIJJJIJIIIIGHFFFFFFEEEEEEEDDDDDFDDDDDDDDD>CDDEDEEDDDDDDCCDDDDDDDDCD    NM:i:0  ms:i:202    AS:i:202    nn:i:0  tp:A:P  cm:i:14 s1:i:94 s2:i:0

--------------------------------

HISEQ:426:C5T65ACXX:5:2302:4953:2090    4   *   0   0   *   *   0   0
AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAAATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAACAACCAAAAAGAGACGAACAA
CCCFFDDFAFFBHJHGGGIHIJBGGHIIJJJJJJJHGEIJGIFIIIHCBGHIJIIIIIJJHHHHEF@D@=;=,0)0&5&))+(((+((((&+(((()&&)(


---------------------------------------

HISEQ:426:C5T65ACXX:5:2302:1943:2127    101 0   101 +   ENA|K01711|K01711.1 15894   447 548 101 101 60  NM:i:0  ms:i:202    AS:i:202    nn:i:0  tp:A:P cm:i:14 s1:i:94 s2:i:0  cg:Z:101M

==> approx-mapping.paf <==
HISEQ:426:C5T65ACXX:5:2302:1943:2127    101 3   97  +   ENA|K01711|K01711.1 15894   450 544 94  94  60  tp:A:P  cm:i:14 s1:i:94 s2:i:0

        """
        pri_only = False
        pattern = "(\d+)([MIDSHNX=])"

        skipped = 0

        with open(self.infile, "r") as fin:
            with open(self.outfile, "w") as fout:
                for lineno, line in enumerate(fin.readlines()):

                    if line.startswith("@"):
                        if line.startswith("@SQ"):
                            #var name = (m = /\tSN:(\S+)/.exec(line)) != null? m[1] : null;
                            #var l = (m = /\tLN:(\d+)/.exec(line)) != null? parseInt(m[1]) :null;
                            #if (name != null && l != null) len[name] = l; 
                            continue
                        continue
                    #else
                    t = line.split()
                    flag = int(t[1])

                    if (t[9] != "*" and t[10] != "*" and len(t[9]) != len(t[10])):
                        raise ValueError("ERROR at line " + str(i) +
                            ":inconsistent SEQ and QUAL lengths - " +
                            str( t[9].length) + " != " + str(t[10].length))

                    if (t[2] == '*' or (flag and 4)):
                        skipped += 1
                        continue

                    if (pri_only and (flag and 0x100)):
                        continue

                    # ??? ???
                    tlen = len(t[2])

                    if (tlen == 0) :
                        raise ValueError("ERROR at line " + str(i) + ": can't find the length of contig " + str(t[2]))
                    match = re.findall("\tnn:i:(\d+)", line)
                    if match:
                        nn = int(match[0])
                    else:
                        nn = 0
                    m = re.findall("\tNM:i:(\d+)", line)
                    if match:
                        NM = int(match[0])
                        have_NM = True
                    else:
                        NM = 0
                        have_NM = False

                    # nn = (m = /\tnn:i:(\d+)/.exec(line)) != null? parseInt(m[1]) : 0;
                    # NM = (m = /\tNM:i:(\d+)/.exec(line)) != null? parseInt(m[1]) : null;

                    # have_NM = NM == null? false : true
                    NM += nn
                    clip = [0, 0]
                    I = [0, 0]
                    D = [0, 0]
                    M = 0
                    N = 0
                    ql, tl, mm = 0, 0, 0,
                    ext_cigar = False
                    n_cigar = 0
                    for count, letter in re.findall(pattern, t[5]):
                        l = int(count)
                        if (letter == 'M'):
                            M += l
                            ql += l
                            tl += l
                            ext_cigar = False
                        elif (letter == 'I'):
                            I[0] += 1
                            I[1] += l
                            ql += l
                        elif (letter == 'D'):
                            D[0] += 1
                            D[1] += l
                            tl += l
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

                    if (n_cigar > 65535):
                         logger.warning("at line " + str(lineno) + ": " + n_cigar + " CIGAR operations")

                    #if (tl + int(t[3]) - 1 > tlen):
                    #    print(tl, t[3], tlen)
                    #    logger.warning("line " + str(lineno) + ": alignment end position larger than ref length; skipped");
                    #    continue

                    if (t[9] != '*' and len(t[9]) != ql) :
                        logger.warning("at line " + str(lineno) + ": SEQ length inconsistent with CIGAR(" + len(t[9]) + " != " + ql + "); skipped")
                        continue

                    if (have_NM is False or ext_cigar):  # to check
                         NM = I[1] + D[1] + mm

                    if (NM < I[1] + D[1] + mm):
                        logger.warning("at line " + str(lineno) + ": NM is less than the total number of gaps (" + str(NM) + " < " + str(I[1]+D[1]+mm) + ")")
                        NM = I[1] + D[1] + mm

                    extra = ["mm:i:"+ str(NM-I[1]-D[1]), "io:i:"+str(I[0]), "in:i:"+str(I[1]),"do:i:"+str(D[0]), "dn:i:"+ str(D[1])]
                    match = M - (NM - I[1] - D[1]);
                    blen = M + I[1] + D[1]
                    qlen = M + I[1] + clip[0] + clip[1]
                    if (flag and 16):
                        qs = clip[1]
                        qe = qlen - clip[0]
                    else:
                        qs = clip[0]
                        qe = qlen - clip[1]
                    ts = int(t[3]) - 1
                    te = ts + M + D[1] + N

                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
                    #FIXME
                    tlen = -1
                    a = [t[0], qlen, qs, qe, "-" if flag and 16 else '+', t[2], tlen, ts, te, match, blen, t[4]]
                    a = [str(x) for x in a]
                    fout.write("\t".join(a) + "\t" + "\t".join(extra) + "\n")

        print(skipped)
