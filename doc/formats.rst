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


PAF
---------

:reference: https://github.com/lh3/miniasm/blob/master/PAF.md

PAF is a text format describing the approximate mapping positions between two
set of sequences. PAF is TAB-delimited with each line consisting of the
following predefined fields:

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
