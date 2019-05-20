.. _formats:


Here below we provide a description of all formats used in bioconvert as well as
other less known formats not included. We annotate the formats:

- Type: sequence, assembly, alignement, other
- Format: binary, human-readable.
- Status: deprecated or not included

..
    .. admonition:: Bioconvert conversions
    .. seealso::  
    .. admonition::  References

Formats
==========

Here below, we provide a list of formats used in bioinformatics or computational
biology. Most of these formats are used in **Bioconvert** and available for
conversion to another formats. Some are available for book-keeping. 

We hope that this page will be useful to all developers and scientists. Would
you like to contribute, please edit the file in our github **doc/formats.rst**.



.. _format_twobit

.2bit (twobit)
--------------

:Format: binary
:Status: available
:Type: sequence


A **2bit** file stores multiple DNA sequences (up to 4 Gb total) in a compact
randomly-accessible format. The file contains masking information as well as the
DNA itself.

The file begins with a 16-byte header containing the following fields:

  - signature: the number 0x1A412743 in the architecture of the machine that created the file
  - version: zero for now. Readers should abort if they see a version number higher than 0
  - sequenceCount: the number of sequences in the file 
  - reserved: always zero for now

All fields are 32 bits unless noted. If the signature value is not as given, the
reader program should byte-swap the signature and check if the swapped version
matches. If so, all multiple-byte entities in the file will have to be
byte-swapped. This enables these binary files to be used unchanged on different
architectures.

The header is followed by a file index, which contains one entry for each
sequence. Each index entry contains three fields:

    - nameSize: a byte containing the length of the name field
    - name: the sequence name itself (in ASCII-compatible byte string), of variable length depending on nameSize
    - offset: the 32-bit offset of the sequence data relative to the start of the file, not aligned to any 4-byte padding boundary

The index is followed by the sequence records, which contain nine fields:

    - dnaSize - number of bases of DNA in the sequence
    - nBlockCount - the number of blocks of Ns in the file (representing unknown sequence)
    - nBlockStarts - an array of length nBlockCount of 32 bit integers indicating the (0-based) starting position of a block of Ns
    - nBlockSizes - an array of length nBlockCount of 32 bit integers indicating the length of a block of Ns
    - maskBlockCount - the number of masked (lower-case) blocks 
    - maskBlockStarts - an array of length maskBlockCount of 32 bit integers indicating the (0-based) starting position of a masked block
    - maskBlockSizes - an array of length maskBlockCount of 32 bit integers indicating the length of a masked block
    - reserved - always zero for now
    - packedDna - the DNA packed to two bits per base, represented as so: T - 00, C - 01, A - 10, G - 11. The first base is in the most significant 2-bit byte; the last base is in the least significant 2 bits. For example, the sequence TCAG is represented as 00011011.

.. admonition:: Reference:

    - http://genome.ucsc.edu/FAQ/FAQformat.html#format7








.. _format_abi:

ABI
---
:Format: binary
:Status: available
:Type: sequence

ABI are trace files that include the PHRED quality scores for the base calls.
This allows ABI to FASTQ conversion. Note that each ABI file contains one and only one sequence (no need for indexing the file). The trace data contains probablities of the four nucleotide bases along the sequencing run together with the sequence deduced from that data. ABI trace is a binary format.

File format produced by ABI sequencing machine. It produces ABI "Sanger" capillary sequence

.. admonition:: Bioconvert conversions:

    :class:`~bioconvert.abi2qual.ABI2QUAL`,
    :class:`~bioconvert.abi2fastq.ABI2FASTQ`,
    :class:`~bioconvert.abi2fasta.ABI2FASTA`

.. seealso:: :ref:`scf`, :class:`~bioconvert.scf2fasta.SCF2Fasta`,
    :class:`~bioconvert.scf2fastq.SCF2Fastq`,

.. admonition::  References

    - ftp://saf.bio.caltech.edu/pub/software/molbio/abitools.zip
    - https://github.com/jkbonfield/io_lib/
    - http://www6.appliedbiosystems.com/support/software_community/ABIF_File_Format.pdf


.. _format_asqg:

ASQG
----

:Format: human-readable
:Status: not included (deprecated)
:Type: assembly

The ASQG format describes an assembly graph. Each line is a tab-delimited
record. The first field in each record describes the record type. The three
types are:

- HT: Header record. This record contains metadata tags for the file version
  (VN tag) and parameters associated with the graph (for example the minimum
  overlap length).
- VT: Vertex records. The second field contains the vertex identifier, the
  third field contains the sequence. Subsequent fields contain optional tags.
- ED: Edge description records. Fields are:
    - sequence 1 name
    - sequence 2 name
    - sequence 1 overlap start (0 based)
    - sequence 1 overlap end (inclusive)
    - sequence 1 length
    - sequence 2 overlap start (0 based)
    - sequence 2 overlap end (inclusive)
    - sequence 2 length
    - sequence 2 orientation (1 for reversed with respect to sequence 1)
    - number of differences in overlap (0 for perfect overlaps, which is the default).

Example::

    HT  VN:i:1  ER:f:0  OL:i:45 IN:Z:reads.fa   CN:i:1  TE:i:0
    VT  read1   GATCGATCTAGCTAGCTAGCTAGCTAGTTAGATGCATGCATGCTAGCTGG
    VT  read2   CGATCTAGCTAGCTAGCTAGCTAGTTAGATGCATGCATGCTAGCTGGATA
    VT  read3   ATCTAGCTAGCTAGCTAGCTAGTTAGATGCATGCATGCTAGCTGGATATT
    ED  read2 read1 0 46 50 3 49 50 0 0
    ED  read3 read2 0 47 50 2 49 50 0 0

.. admonition:: References

    - https://github.com/jts/sga/wiki/ASQG-Format


.. _format_bai:

BAI
---
:Format: binary
:Status: not included
:Type: index

The index file of a BAM file is a BAI file format. The BAI files are 
not used in **Bioconvert**. 


.. _format_bam:

BAM
---

:Format: binary
:Status: included
:Type: Sequence alignement

The BAM (Binary Alignment Map) is the binary version of the Sequence 
Alignment Map (:ref:`format_sam`) format. It is a compact and index-able representation 
of nucleotide sequence alignments. 

.. admonition:: Bioconvert Conversions

    :class:`~bioconvert.bam2sam.BAM2SAM`,
    :class:`~bioconvert.bam2cram.BAM2CRAM`,
    :class:`~bioconvert.bam2bedgraph.BAM2BEDGRAPH`,
    :class:`~bioconvert.bam2bed.BAM2BED`,
    :class:`~bioconvert.bam2bigwig.BAM2BIGWIG`,
    :class:`~bioconvert.bam2fasta.BAM2FASTA`,
    :class:`~bioconvert.bam2fastq.BAM2FASTQ`,
    :class:`~bioconvert.bam2json.BAM2JSON`,
    :class:`~bioconvert.bam2tsv.BAM2TSV`,
    :class:`~bioconvert.bam2wiggle.BAM2WIGGLE`

.. admonition:: References

    - http://samtools.github.io/hts-specs/SAMv1.pdf
    - http://genome.ucsc.edu/goldenPath/help/bam.html

.. seealso:: The :ref:`format_sam` and :ref:`format_bai` formats.


.. _format_bcf:

BCF
---

:Format: binary
:Status: included
:Type: variant

Binary version of the Variant Call Format (VCF).

.. admonition:: Bioconvert conversions

    - :class:`~bioconvert.bcf2vcf.BCF2VCF`
    - :class:`~bioconvert.vcf2bcf.VCF2BCF`





.. _format_bcl:

BCL
---

:Format: binary
:Status: not included
:Type: sequence

BCL is the raw format used by Illumina sequencer. This data is converted into
:ref:`FastQ  <format_fastq>` thanks to a tool called bcl2fastq. This type of conversion is not included
in **Bioconvert**


.. _format_bedgraph:

BEDGRAPH
--------

.. _format_bed:

BED
---

:reference: http://genome.ucsc.edu/FAQ/FAQformat.html#format1

BED file must has at least 3 columns (chrom, start, end).

.. _format_bigwig:

BIGWIG
------


.. _format_cram:

CRAM
----

:Format: binary
:Status: not included
:Type: Alignment

The CRAM file format is a more dense form of BAM files with the benefit of
saving much disk space. While BAM files contain all sequence data within a file,
CRAM files are smaller by taking advantage of an additional external reference
sequence file. This file is needed to both compress and decompress the read
information.

.. seealso:: :ref:`format_bam`


.. admonition:: Bioconvert Conversions

    - :class:`~bioconvert.bam2sam.BAM2CRAM`
    - :class:`~bioconvert.bam2cram.SAM2CRAM`
    - :class:`~bioconvert.bam2sam.CRAM2BAM`
    - :class:`~bioconvert.bam2cram.CRAM2SAM`


.. _format_csv:

CSV
---

:Type: database
:Status: included

A comma-separated values format is a delimited text file that uses a
comma to separate values. See :ref:`format_csv` format page for
details.

.. admonition:: References

    - https://en.wikipedia.org/wiki/Comma-separated_values



.. _format_fasta:

FastA
-----

:Format: human-readable
:Status: included
:Type: Sequence

This refers to the input FASTA file format where each record starts
with a ">" line. Resulting sequences have a generic alphabet by default. 
There is no standard file extension for a text file containing FASTA formatted sequences. Although
their is a plethora of ad-hoc file extensions: fasta, fas, fa, seq, fsa, fna, ffn, faa, frn, we use only fasta, fa and fst within **Bioconvert**.


.. admonition:: Bioconvert conversions

    - :class:`~bioconvert.fastq2fasta.FastQ2FastA`
    - :class:`~bioconvert.fasta2fasta.FastA2FastQ`
    - :class:`~bioconvert.fasta2clustal.FastA2Clustal`
    - :class:`~bioconvert.fasta2nexus.FastA2Nexus`
    - :class:`~bioconvert.fasta2twobit.FastA2TwoBit`

.. seealso:: :ref:`format_fastq` and :ref:`format_qual`
.. admonition::  References

    -  http://en.wikipedia.org/wiki/FASTA_format


.. _format_fastg:

FastG
-----

:Format:
:Status: not included 
:Type: assembly


:reference: http://fastg.sourceforge.net/FASTG_Spec_v1.00.pdf


.. _format_fastq:

FastQ
-----

:Format: human-readable
:Status: included
:Type: Sequence

FASTQ files include sequences in :ref:`format_fasta` format and their 
qualities (:ref:`format_qual`). In general, *fastq*
refers to Sanger style FASTQ files which encode PHRED qualities using an
ASCII offset of 33. See also the incompatible "fastq-solexa" and "fastq-illumina"
variants used in early Solexa/Illumina pipelines, Illumina pipeline 1.8 produces Sanger FASTQ.
Be aware that there are different FASTQ formats for different sequencing technologiess

.. admonition:: Bioconvert conversions

    - :class:`~bioconvert.fastq2fasta.FastQ2FastA`
    - :class:`~bioconvert.fasta2fasta.FastA2FastQ`

.. seealso:: :ref:`format_fasta` and ref:`format_qual`

.. _gfa_format:

GFA
---

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

.. _format_json:

JSON
----

TODO

.. _format_nexus:

Nexus
-----------

The NEXUS multiple alignment format, also known as PAUP format. 



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
- father ID: 0 if unknown. If specified, must also appear as an individual in the file
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


.. _format_qual:

QUAL
----

:Format: human-readable
:Status: included
:Type: Sequence

QUAL files include qualities of each nucleotide in :ref:`format_fasta` format.

.. admonition:: Bioconvert conversions

    - :class:`~bioconvert.fastq2fasta.FastQ2FastA`
    - :class:`~bioconvert.fasta2fasta.FastA2FastQ`

.. seealso:: :ref:`format_fasta` and :ref:`format_fastq`


BED for plink
~~~~~~~~~~~~~~
This BED format  is the binary PED file. Not to be confused with BED format used
with BAM files.

BIM files
~~~~~~~~~

The fields are 

- chromosome number (integer)
- SNP marker ID (string)
- SNP generit position (cM) (float)
- SNP physical position (bp)
- Allele 1
- Allele 2

So, it is like the MAP with the 2 alleles, and the format is binary.

.. _format_fam:

FAM 
~~~

The first 6 columns of the PED file.


.. _format_sam:

SAM
---

:Format: human readable
:Status: included
:Type: alignment


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

The optional fields are tool-dependent. For instance with BWA mapper, we can get these tags

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




.. admonition::  References

    - http://samtools.github.io/hts-specs/SAMv1.pdf
    - http://genome.ucsc.edu/goldenPath/help/bam.html


.. _format_scf:

Trace File Format - Sequence Chromatogram Format (SCF)
------------------------------------------------------

:reference: https://wiki.nci.nih.gov/display/TCGA/Sequence+trace+files
:reference: http://staden.sourceforge.net/manual/formats_unix_2.html

Trace files are binary files containing raw data output from automated sequencing instruments.
This converter was converted from BioPerl.


SCF file organisation (more or less)

====================================== ====================================
Length in bytes                        Data
====================================== ====================================
128                                    header
Number of samples * sample size        Samples for A trace
Number of samples * sample size        Samples for C trace
Number of samples * sample size        Samples for G trace
Number of samples * sample size        Samples for T trace
Number of bases * 4                    Offset into peak index for each base
Number of bases                        Accuracy estimate bases being 'A'
Number of bases                        Accuracy estimate bases being 'C'
Number of bases                        Accuracy estimate bases being 'G'
Number of bases                        Accuracy estimate bases being 'T'
Number of bases                        The called bases
Number of bases * 3                    Reserved for future use
Comments size                          Comments
Private data size                      Private data
====================================== ====================================

.. _format_tsv: 

TSV
---

:Type: database
:Status: included

A tab-separated values format is a delimited text file that uses a
tab character to separate values. See :ref:`format_csv` format page for
details.


.. admonition:: Bioconvert conversions:

    :class:`~bioconvert.tsv2csv.TSV2CSV`,



.. admonition:: References

    - https://en.wikipedia.org/wiki/Comma-separated_values






Stockholm
---------

The Stockholm alignment format is also known as PFAM format.   


.. _format_vcf:

VCF
---

:Format: human readable
:Status: included
:Type: variant


Variant Call Format (VCF) is a flexible and extendable format for 
storing variation in sequences such as single nucleotide variants,
insertions/deletions, copy number variants and structural variants. 

.. admonition:: Bioconvert conversions:

    - :class:`~bioconvert.bcf2vcf`
    - :class:`~bioconvert.bcf2wiggle`
    - :class:`~bioconvert.vcf2bcf`
    - :class:`~bioconvert.vcf2bed`
    - :class:`~bioconvert.vcf2wiggle`
    - :class:`~bioconvert.vcf2plink`
    - :class:`~bioconvert.vcf2bplink`




Wiggle Track format (WIG)
-------------------------

:reference: http://genome.ucsc.edu/goldenPath/help/wiggle.html

The bigWig format is used for graphing track needs. The wiggle (WIG) format is
an older format for display of dense, continuous data such as GC percent. 
Wiggle data elements must be equally sized. 

Similar format such as the bedGraph format is also an older format used to display sparse data
or data that contains elements of varying size.

For speed and efficiency, wiggle data is compressed with a minor loss of precision when
data is exported from a wiggle track.

.. _format_xls:

XLS
---

:Type: database
:Status: included

Spreadsheet file format (Microsoft Excel file format). 

Until 2007, Microsoft Excel used a proprietary binary file format
called Excel Binary File Format (.XLS). In Excel 2007, the Office Open XML was
introduced. We support the later formnat only.

With bioconvert you can convert an :ref:`format_xls` file into :ref:`format_csv` or :ref:`format_tsv` format. If several
sheets are to be found, you can select one or the other.


.. admonition:: Bioconvert conversions:

    :class:`~bioconvert.xls2csv`,
    :class:`~bioconvert.xlsx2csv`,

.. admonition::  References

    - https://en.wikipedia.org/wiki/Microsoft_Excel#File_formats

.. _format_xlsx:

XLSX
----

:Type: database
:Status: included

Spreadsheet file format in Office Open XML format.


With bioconvert you can convert an :ref:`format_xlsx` file into :ref:`format_csv` or :ref:`format_tsv` format. If several
sheets are to be found, you can select one or the other.


.. admonition:: Bioconvert conversions:

    :class:`~bioconvert.xls2csv`,
    :class:`~bioconvert.xlsx2csv`,

.. seealso::  :ref:`format_xls` format.

.. admonition::  References

    - https://en.wikipedia.org/wiki/Office_Open_XML




TODO
-------
bcf2vcf.py
bcf2wiggle.py
bigbed2wiggle.py
bigwig2bedgraph.py
bigwig2wiggle.py
bplink2plink.py
clustal2fasta.py
clustal2nexus.py
clustal2phylip.py
clustal2stockholm.py
dsrc2gz.py
embl2fasta.py
embl2genbank.py
fasta2clustal.py
fasta2genbank.py
fasta2nexus.py
fasta2phylip.py
fasta2twobit.py
genbank2embl.py
genbank2fasta.py
genbank2gff3.py
gfa2fasta.py
gff22gff3.py
gff3gff2.py
gz2bz2.py
gz2dsrc.py
json2yaml.py
maf2sam.py
newick2nexus.py
newick2phyloxml.py
nexus2clustal.py
nexus2newick.py
nexus2phylip.py
nexus2phyloxml.py
ods2csv.py
phylip2clustal.py
phylip2fasta.py
phylip2nexus.py
phylip2stockholm.py
phylip2xmfa.py
phyloxml2newick.py
phyloxml2nexus.py
plink2bplink.py
plink2vcf.py
sam2paf.py
scf2fasta.py
scf2fastq.py
sra2fastq.py
stockholm2clustal.py
stockholm2phylip.py
tsv2csv.py
twobit2fasta.py
vcf2bcf.py
vcf2bed.py
vcf2bplink.py
vcf2plink.py
vcf2wiggle.py
wig2bed.py
xmfa2phylip.py
yaml2json.py
