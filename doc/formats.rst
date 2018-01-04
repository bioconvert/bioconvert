Formats
==========

ASQG format
--------------

:type: assembly
:reference: https://github.com/jts/sga/wiki/ASQG-Format

The ASQG format describes an assembly graph. Each line is a tab-delimited
irecord. The first field in each record describes the record type. The three
types are:

- HT: Header record. This record contains metadata tags for the file version
(VN tag) and parameters associated with the graph (for example the minimum
overlap length).
- VT: Vertex records. The second field contains the vertex identifier, the
third field contains the sequence. Subsequent fields contain optional tags.
- ED: Edge description records. The second field describes a pair of
overlapping sequences. A full description of this field is below. Subsequent
fields contain optional tags.

Example
~~~~~~~~~~~~~~~~

::

    HT  VN:i:1  ER:f:0  OL:i:45 IN:Z:reads.fa   CN:i:1  TE:i:0
    VT  read1   GATCGATCTAGCTAGCTAGCTAGCTAGTTAGATGCATGCATGCTAGCTGG
    VT  read2   CGATCTAGCTAGCTAGCTAGCTAGTTAGATGCATGCATGCTAGCTGGATA
    VT  read3   ATCTAGCTAGCTAGCTAGCTAGTTAGATGCATGCATGCTAGCTGGATATT
    ED  read2 read1 0 46 50 3 49 50 0 0
    ED  read3 read2 0 47 50 2 49 50 0 0



BAM format
---------------

:reference: http://samtools.github.io/hts-specs/SAMv1.pdf



BED format
---------------

:reference: http://genome.ucsc.edu/FAQ/FAQformat.html#format1

BED file must has at least 3 columns (chrom, start, end).


FastG
----------

:type: assembly
:reference: http://fastg.sourceforge.net/FASTG_Spec_v1.00.pdf




GFA format
-------------

:type: assembly graph
:references: http://gfa-spec.github.io/GFA-spec/,

Overview
~~~~~~~~~~

The Graphical Fragment Assembly (GFA) can be used to represent genome
assemblies. GFA stores sequence graphs as the product of an
assembly, a representation of variation in genomes, splice graphs in genes, or
even overlap between reads from long-read sequencing technology.

The GFA format is a tab-delimited text format for describing a set of sequences
and their overlap. The first field of the line identifies the type of the line.
**Header** lines start with H. **Segment** lines start with S. **Link** lines start with L.
A **containment** line starts with C. A **path** line starts with P.


Terminology
~~~~~~~~~~~~~
- Segment a continuous sequence or subsequence.
- Link an overlap between two segments. Each link is from the end of one segment to the beginning of another segment. The link stores the orientation of each segment and the amount of basepairs overlapping.
- Containment an overlap between two segments where one is contained in the other.
- Path an ordered list of oriented segments, where each consecutive pair of oriented segments are supported by a link record.

See details in the reference above.

Example:
~~~~~~~~~

::

    H   VN:Z:1.0
    S   11  ACCTT
    S   12  TCAAGG
    S   13  CTTGATT
    L   11  +   12  -   4M
    L   12  -   13  +   5M
    L   11  +   13  +   3M
    P   14  11+,12-,13+ 4M,5M


Notes: sometimes you would have extra field (fourth one) on **segment** lines. 
Convertion to fasta will store this fourth line after the name.


GFA version 2
~~~~~~~~~~~~~~~~~~~~~~~~

GFA2 is a generalization of GFA that allows one to specify an assembly graph in
either less detail, e.g. just the topology of the graph, or more detail, e.g.
the multi-alignment of reads giving rise to each sequence. It is further
designed to be a able to represent a string graph at any stage of assembly, from
the graph of all overlaps, to a final resolved assembly of contig paths with
multi-alignments. Apart from meeting these needs, the extensions also supports
other assembly and variation graph types.

Like GFA, GFA2 is tab-delimited in that every lexical token is separated from
the next by a single tab.


PAF (Pairwise mApping Format)
--------------------------------

:reference: https://github.com/lh3/miniasm/blob/master/PAF.md

PAF is a text format describing the approximate mapping positions between two
set of sequences. PAF is used for instance in **miniasm** tool (see reference
above), an ultrafast de novo assembly for long noisy reads. PAF is TAB-delimited 
with each line consisting of the following predefined fields:

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

If PAF is generated from an alignment, column 10 equals the number of sequence
matches, and column 11 equals the total number of sequence matches, mismatches,
insertions and deletions in the alignment. If alignment is not available, column
10 and 11 are still required but can be approximate.

A PAF file may optionally contain SAM-like typed key-value pairs at the end of
each line.

PLINK flat files (MAP/PED)
-------------------------------

PLINK is a used application for analyzing genotypic data. It can be considered  the de-facto standard of the field. The MAP files describes the SNPs and contains those fields:

- chromosome number (integer)
- SNP marker ID (string)
- SNP generit position (cM) (float)
- SNP physical position (bp)

So it contains L lines with 4 columns. All SNPs must be ordered by physical
position. Example::

    X rs3883674 0 32380
    X rs12218882 0 48172
    9 rs10904045 0 48426
    9 rs10751931 0 49949

The PED (pedigree) file describes the individuals and the genetic data. The PED
file can be spaced or tab delimited. Each line corresponds to a single
individual. The first 6 columns are:

- family ID (or pedigree name): a unique alpha numeric identifier 
- individual ID: should be unique within his family
- father ID: 0 if unknown. If specified, must also appear as an individual in
the file
- mother ID: same as above
- Sex: 1 Male, 2 Female
- Phenotype

- columns 7 and 8 code for the observed alleles at SNP1
- comumns 9 and 10 code for the observed alleles at SNP2 and so on

missing data are coded as "0 0". So we havez N lines 2L + 6 columns where N is
the number of individuals and L the numbers of SNPs

PLINK binary files (BED/BIM/FAM)
-------------------------------------
Same information as plink flat files. 

BED for plink
~~~~~~~~~~~~~~
This BED format  is the binary PED file. Not to be confused with BED format used
with BAM files.

BIM files
~~~~~~~~~~~~~~~~~~~~

The fields are 

- chromosome number (integer)
- SNP marker ID (string)
- SNP generit position (cM) (float)
- SNP physical position (bp)
- Allele 1
- Allele 2

So, it is like the MAP with the 2 alleles, and the format is binary.

FAM files
~~~~~~~~~~~~~~~~~~~~~~~

The first 6 columns of the PED file.






SAM format
-------------

:reference: https://samtools.github.io/hts-specs/SAMv1.pdf


In the SAM format, each alignment line typically represents the linear alignment
of a segment.  Each line has 11 mandatory  fields in the same order. Their values
can be `0` or `*` if the field is unavailable. Here is an overview of those
fields:

======= ======= ======= ======================= ======================================
Col     Field   Type    Regexp/Range            Brief description
======= ======= ======= ======================= ======================================
1       QNAME   String  [!-?A-~]{1,254}         Query template NAME
2       FLAG    Int     [0,2^16-1]              bitwise FLAG
3       RNAME   String  \*|[!-()+-<>-~][!-~]*   Reference sequence NAME
4       POS     Int     [0,2^31-1]              1-based leftmost mapping POSition
5       MAPQ    Int     [0,2^8-1]               MAPping Quality
6       CIGAR   String  \*|([0-9]+[MIDNSHPX=])+ CIGAR string
7       RNEXT   String  \*|=|[!-()+-<>-~][!-~]* Ref.  name of the mate/next read
8       PNEXT   Int     [0,2^31-1]              Position of the mate/next read
9       TLEN    Int     [-2^31+1,2^31-1]        observed Template LENgth
10      SEQ     String  \*|[A-Za-z=.]+          segment SEQuence
11      QUAL    String  [!-~]+                  ASCII of Phred-scaled base QUALity+33
======= ======= ======= ======================= ======================================

All  optional   fields  follow  the TAG:TYPE:VALUE format  where TAG is  a  two-character  string  that  matches /[A-Za-z][A-Za-z0-9]/ .  Each TAG can only appear once in one alignment line.

The tag `NM:i:2` means: Edit distance to the reference (number of changes
necessary to make this equal to the reference, exceluding clipping).


The optional fields are tool-dependent. 

From BWA documentation, we can get this

==== ==================================================
Tag         Meaning
==== ==================================================
NM         Edit distance
MD         Mismatching positions/bases
AS         Alignment score
BC         Barcode sequence
X0         Number of best hits
X1         Number of suboptimal hits found by BWA
XN         Number of ambiguous bases in the referenece
XM         Number of mismatches in the alignment
XO         Number of gap opens
XG         Number of gap extentions
XT         Type: Unique/Repeat/N/Mate-sw
XA         Alternative hits; format: (chr,pos,CIGAR,NM;)*
XS         Suboptimal alignment score
XF         Support from forward/reverse alignment
XE         Number of supporting seeds
==== ==================================================

Note that XO and XG are generated by BWT search while the CIGAR string by
Smith-Waterman alignment. These two tags may be inconsistent with the CIGAR
string. This is not a bug

`SA:Z`: Other canonical alignments in a chimeric alignment, in the format of: (rname,pos,strand,CIGAR,mapQ,NM;)+. Each element in the semi-colon delimited list represents a part of the chimeric alignment. Conventionally, at a supplementary line, the first element points to the primary line.



Wiggle Track format (WIG)
------------------------------

:reference: http://genome.ucsc.edu/goldenPath/help/wiggle.html

The bigWig format is used for graphing track needs. The wiggle (WIG) format is
an older format for display of dense, continuous data such as GC percent. 
Wiggle data elements must be equally sized. 

Similar format such as the bedGraph format is also an older format used to display sparse data
or data that contains elements of varying size.

For speed and efficiency, wiggle data is compressed with a minor loss of precision when
data is exported from a wiggle track.
