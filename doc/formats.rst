

.. _formats:

Formats
=======

Here below, we provide a list of formats used in bioinformatics or computational
biology. Most of these formats are used in **Bioconvert** and available for
conversion to another formats. Some are available for book-keeping.

We hope that this page will be useful to all developers and scientists. Would
you like to contribute, please edit the file in our github **doc/formats.rst**.

If you wish to update this page, please see the :ref:`developer_guide` page.


.. - Type: sequence, assembly, alignement, other, index, variant, database,
   compression
.. - Format: binary, human-readable
.. - Status: deprecated, included, not included



.. _format_twobit:

TWOBIT
------

:Format: binary
:Status: available
:Type: sequence


A **2bit** file stores multiple DNA sequences (up to 4 Gb total) in a compact
randomly-accessible format. The file contains masking information as well as the
DNA itself.

The file begins with a 16-byte header containing the following fields:

  - signature: the number 0x1A412743 in the architecture of the machine that
    created the file
  - version: zero for now. Readers should abort if they see a version number
    higher than 0
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
    - name: the sequence name itself (in ASCII-compatible byte string), of
      variable length depending on nameSize
    - offset: the 32-bit offset of the sequence data relative to the start of
      the file, not aligned to any 4-byte padding boundary

The index is followed by the sequence records, which contain nine fields:

    - dnaSize - number of bases of DNA in the sequence
    - nBlockCount - the number of blocks of Ns in the file (representing
      unknown sequence)
    - nBlockStarts - an array of length nBlockCount of 32 bit integers
      indicating the (0-based) starting position of a block of Ns
    - nBlockSizes - an array of length nBlockCount of 32 bit integers
      indicating the length of a block of Ns
    - maskBlockCount - the number of masked (lower-case) blocks
    - maskBlockStarts - an array of length maskBlockCount of 32 bit integers
      indicating the (0-based) starting position of a masked block
    - maskBlockSizes - an array of length maskBlockCount of 32 bit integers
      indicating the length of a masked block
    - reserved - always zero for now
    - packedDna - the DNA packed to two bits per base, represented as
      so: T - 00, C - 01, A - 10, G - 11. The first base is in the most
      significant 2-bit byte; the last base is in the least significant 2 bits.
      For example, the sequence TCAG is represented as 00011011.

.. admonition:: Bioconvert conversions

    :class:`~bioconvert.twobit2fasta.TWOBIT2FASTA`

.. admonition:: Reference:

    - http://genome.ucsc.edu/FAQ/FAQformat.html#format7


.. _format_agp:

AGP
---
:Format: human-readable
:Status: 
:Type: assembly

AGP files are used to describe the assembly of a sequences from smaller
fragments. The large object can be a contig, a scaffold (supercontig), or a chromosome. Each
line (row) of the AGP file describes a different piece of the object, and has
the column entries defined below. Several format exists: 1.0, 2.0, 2.1


you can validate your AGP file using this website:
https://www.ncbi.nlm.nih.gov/projects/genome/assembly/agp/agp_validate.cgi

.. admonition::  References

    - https://www.ebi.ac.uk/ena/submit/agp-files
    - https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Validation/


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

.. seealso:: :ref:`format_scf`, :class:`~bioconvert.scf2fasta.SCF2FASTA`,
    :class:`~bioconvert.scf2fastq.SCF2FASTQ`,

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
    :class:`~bioconvert.bam2bed.BAM2COV`,
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

Binary version of the Variant Call Format (:ref:`format_vcf`).

.. admonition:: Bioconvert conversions

    :class:`~bioconvert.bcf2vcf.BCF2VCF`, :class:`~bioconvert.vcf2bcf.VCF2BCF`.
    :class:`~bioconvert.bcf2wiggle.BCF2WIGGLE`


.. _format_bcl:

BCL
---

:Format: binary
:Status: not included
:Type: sequence

BCL is the raw format used by Illumina sequencer. This data is converted into
:ref:`FastQ  <format_fastq>` thanks to a tool called bcl2fastq. This type of conversion is not included
in **Bioconvert**. Indeed,  Illumina provides a **bcl2fastq** executable and its user guide is available online. In most cases, the BCL files are already converted and users will only get the FastQ files so we will not provide such converter.

.. admonition:: References

    - https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq_letterbooklet_15038058brpmi.pdf
    - http://bioinformatics.cvr.ac.uk/blog/how-to-demultiplex-illumina-data-and-generate-fastq-files-using-bcl2fastq/


BED for plink
-------------

:Format: binary
:Status: included
:Type: genotypic

This BED format  is the binary PED file. Not to be confused with BED format used
with BAM files. Please see :ref:`format_plink_binary` section.


.. _format_bedgraph:

BEDGRAPH
--------

:Format: human-readable
:Status: included
:Type: database

BedGraph is a subset of BED12 format. It is a 4-columns tab-delimited file with
chromosome name, start and end positions and the fourth column is a number that is
often used to show coverage depth. So, this is the same format as the
:ref:`format_bed4` format.
Example::

    chr1    0     75  0
    chr1    75   176  1
    chr1    176  177  2


.. seealso:: :ref:`format_bed`

.. _format_bed:

BED
---

:Format: human-readable
:Status: not included
:Type: database

A Browser Extensible Data (BED) file is a tab-delimited text file. It is a
concise way to represent genomic features and annotations.

The BED file is a very versatile format, which makes it difficult to handle in **Bioconvert**. So, let us describe exhaustively the BED format.

Although the BED description format supports up to 12 columsn, only the first 3
are required for some tools such as the UCSC browser, Galaxy, or bedtools
software.

So, in general BED lines have 3 required fields and nine additional
optional fields.

Generally, all BED files have the same extensions (.bed) irrespective of the
number of columns. We can refer to the 3-columns version as BED3, the 4-columns 
BED as :ref:`format_bed4` and so on.

The number of fields per line must be consistent. If some fields are empty,
additional column information must be filled for consistency (e.g., with a ".").
BED fields can be whitespace-delimited or tab-delimited although some
variations of BED types such as "bed Detail" require a tab character
delimitation for the detail columns (see Note box here below).


.. note:: *BED detail* format

    It is an extension of BED format plus 2 additional fields.
    The first one is an ID, which can be used in place of the name field
    for creating links from the details pages. The second additional field
    is a description of the item, which can be a long description and can
    consist of html.

    Requirements:
        - fields must be tab-separated
        - "type=bedDetail" must be included in the track line,
        - the name and position fields should uniquely describe items
          so that the correct ID and description will be displayed on
          the details pages.

     The following example uses the first 4 columns of BED format,
     but up to 12 may be used. Note the header, which contains the
     type=bedDetail string.::

         track name=HbVar type=bedDetail description="HbVar custom track" db=hg19  visibility=3 url="blabla.html"
         chr11  5246919 5246920 Hb_North_York   2619    Hemoglobin variant
         chr11  5255660 5255661 HBD c.1 G>A 2659    delta0 thalassemia
         chr11  5247945 5247946 Hb Sheffield    2672    Hemoglobin variant
         chr11  5255415 5255416 Hb A2-Lyon  2676    Hemoglobin variant
         chr11  5248234 5248235 Hb Aix-les-Bains    2677    Hemoglobin variant


.. warning:: Browser such as the Genome Browser (http://genome.ucsc.edu/) can visualise BED
    files. Usually, BED files can be annotated using header lines, which begin with the
    word "browser" or "track" to assist the browser in the display and interpretation.

    Such annotation track header lines are not permissible in utilities such as
    bedToBigBed, which convert lines of BED text to indexed binary files.


The file description below is modified from: http://genome.ucsc.edu/FAQ/FAQformat#format1.

The first three required BED fields are:

1. **chrom** - The name of the chromosome (e.g. chr3) or scaffold.
2. **chromStart** - The starting position of the feature in the chromosome.
   The first base in a chromosome is numbered 0.
3. **chromEnd** - The ending position of the feature in the chromosome.
   The chromEnd base is not included in the display of the feature.

The 9 additional optional BED fields are:

4. **name** - Label of the BED line

5. **score** - A score between 0 and 1000. In Genome Browser, the track line
   useScore attribute is set to 1 for this annotation data set, the score value
   will determine the level of gray in which this feature is displayed.

6. **strand** - Defines the strand. Either "." (=no strand) or "+" or "-".

7. **thickStart** - The starting position at which the feature is drawn thickly
   (for example, the start codon in gene displays). When there is no thick part,
   thickStart and thickEnd are usually set to the chromStart position.

8. **thickEnd** - The ending position at which the feature is drawn thickly.

9. **itemRgb** - An RGB value of the form R,G,B (e.g. 255,0,0).

10. **blockCount** - The number of blocks (exons) in the BED line.

11. **blockSizes** - A comma-separated list of the block sizes.
    The number of items in this list should correspond to blockCount.

12. **blockStarts** - A comma-separated list of block starts. Should be
    calculated relative to chromStart. The number of items in this list
    should correspond to blockCount.

In BED files with block definitions, the first blockStart value must be 0, so that the first block begins at chromStart. Similarly, the final blockStart position plus the final blockSize value must equal to chromEnd. Blocks may not overlap.

Here is a simple example::

    track name=pairedReads description="Clone Paired Reads" useScore=1
    chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512
    chr22 2000 6000 cloneB 900 - 2000 6000 0 2 433,399, 0,3601


.. note:: If your data set is BED-like, but it is very large (over 50MB)
    you can convert it to a :ref:`format_bigbed` format.

.. seealso:: :ref:`format_bedgraph`


.. _format_bed3:

BED3
----

A BED3 is supported by bedtools. It is a BED file where each feature is
described by chrom, start and end (with tab-delimited values). Example::

    chr1    100    120

See :ref:`format_bed` section for details.

.. _format_bed4:

BED4
----

A BED4 is a :ref:`format_bed` file where each feature is
described by chrom, start, end and name (with tab-delimited values). The last 
column could also be a number. Example::

    chr1    100    120    gene1

See :ref:`format_bed` section for details.

.. seealso:: :ref:`format_bedgraph`

.. _format_bed5:

BED5
----

A BED5 is supported by bedtools. It is a BED file where each feature is
described by chrom, start, end, name and score(with tab-delimited values). Example::

    chr1    100    120    gene1 0

See :ref:`format_bed` section for details.


.. _format_bed6:

BED6
----

A BED6 is supported by bedtools. It is a BED file where each feature is
described by chrom, start, end, name, score and strand (with tab-delimited values). Example::

    chr1    100    120    gene1 0 +

See :ref:`format_bed` section for details.


.. _format_bed12:

BED12
-----

A BED12 is supported by bedtools. It is a BED file where each feature is
described by all 12 BED fields. Example::

    chr1    100    120    gene1 0 + 100 100 0 3 1,2,3 4,5,6

See :ref:`format_bed` section.


.. _format_bigbed:

BIGBED
------

:Format: binary
:Status: included
:Type: database/track


The **bigBed** format stores annotation items. BigBed files are created initially from BED type files. The resulting bigBed files are in an indexed binary format. The main advantage of the bigBed files is that only the portions of the files needed to display a particular region is used.

.. admonition:: bioconvert conversions

    :class:`~bioconvert.bigbed2bed.BIGBED2COV`, 
    :class:`~bioconvert.bigbed2wiggle.BIGBED2WIGGLE`

.. admonition:: References

    - http://genome.ucsc.edu/goldenPath/help/bigBed.html
    - https://github.com/deeptools/pyBigWig


.. _format_bigwig:

BIGWIG
------

:Format: binary
:Status: included
:Type: database/track

The bigWig format is useful for dense, continuous data. They can be created from
wiggle file (:ref:`format_wiggle`). This type of file is an indexed binary format.

Wiggle data must be continuous unlike :ref:`format_bed`. You can convert a
BED/BEDGraph to bigwig using :class:`~bioconvert.bedgraph2bigwig.BEDGRAPH2BIGWIG`.


To create a bigwig from a wiggle, yo need to remove the existing "track" header

.. admonition:: Bioconvert conversions::

    :class:`~bioconvert.bigwig2wiggle.BIGWIG2WIGGLE`,
    :class:`~bioconvert.bigwig2wiggle.bedgraph2bigwig.BEDGRAPH2BIGWIG`




.. note:: Wiggle, bigWig, and bigBed files use 0-based half-open coordinates, which are
    also used by this extension. So to access the value for the first base on chr1,
    one would specify the starting position as 0 and the end position as 1.
    Similarly, bases 100 to 115 would have a start of 99 and an end of 115. This is
    simply for the sake of consistency with the underlying bigWig file and may
    change in the future in various formats and tools dealing with those formats.


.. admonition:: References:

    - https://genome.ucsc.edu/goldenpath/help/bigWig.html


.. _format_bim:

BIM
---

:Format: human-readable
:Status: included
:Type: variants

The BIM formatted file is a variant information file accompanying
a .bed or biallelic .pgen binary genotype table. Please see :ref:`format_plink_binary` section.


The fields are:

- chromosome number (integer)
- SNP marker ID (string) / variant ID
- SNP generic position (cM) (float) / position in centimorgans (safe to use dummy value 0)
- SNP physical position (bp) (1-based)
- Alternate allele code
- Reference allele code

Here is an example::

       1    rs0     0   1000    0   1
       1    rs10    0   1001    2   1


.. _format_bz2:

BZ2
---

:Format: binary
:Status: included
:Type: Compression


**bzip2** is a file compression program that uses the Burrows–Wheeler algorithm. Extension is usually .bz2
The BZ2 compression is usually better than gzip for Fastq format compression (factor 2-3).

.. admonition:: Bioconvert conversions:

    :class:`~bioconvert.gz2bz2`,
    :class:`~bioconvert.gz2dsrc`
    :class:`~bioconvert.bz22gz`,
    :class:`~bioconvert.dsrc2gz`

.. _format_cov:

COV
---

A simple TSV file with 3 columns to store coverage in a continuous way. First
column is contig/chromosome name, second is position and third is coverage.
Expected positions are continuous. The :ref:`format_bedgraph` stores an extra
column but can be a more compact way of storing coverage/depth.

Example::

    chr1   1    10
    chr1   2    11
    chr1   3    15
    chr1   4    12
    chr1   5    11


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

    :class:`~bioconvert.bam2sam.BAM2CRAM`, :class:`~bioconvert.bam2cram.SAM2CRAM`,
    :class:`~bioconvert.bam2sam.CRAM2BAM`, :class:`~bioconvert.bam2cram.CRAM2SAM`.


.. _format_clustal:

CLUSTAL
-------

:Format: human-readable
:Status: included
:Type: multiple alignment


In a Clustal format, the first line in the file must start with the words "CLUSTAL W"
or "CLUSTALW". Nevertheless, many such files starts with CLUSTAL or CLUSTAL X.
Other information in the first line is ignored. One or more empty lines.
One or more blocks of sequence data. Each block consists of one line for each sequence
in the alignment. Each line consists of the sequence name white space up to 60 sequence symbols.
optional - white space followed by a cumulative count of residues for the  sequences
A line showing the degree of conservation for the columns of the alignment in this block.
One or more empty lines.

Some rules about representing sequences:

- Case does not matter.
- Sequence symbols should be from a valid alphabet.
- Gaps are represented using hyphens ("-").
- The characters used to represent the degree of conservation are
    - `*`  - : all residues or nucleotides in that column are identical
    - `:`  - : conserved substitutions have been observed
    - `.`  - : semi-conserved substitutions have been observed
    - <SPACE>  - : no match.

Here is an example of a multiple alignment in CLUSTAL W format::

    CLUSTAL W (1.82) multiple sequence alignment


    FOSB_MOUSE      MFQAFPGDYDSGSRCSSSPSAESQYLSSVDSFGSPPTAAASQECAGLGEMPGSFVPTVTA 60
    FOSB_HUMAN      MFQAFPGDYDSGSRCSSSPSAESQYLSSVDSFGSPPTAAASQECAGLGEMPGSFVPTVTA 60
                    ************************************************************

    FOSB_MOUSE      TSSFVLTCPEVSAFAGAQRTSGSEQPSDPLNSPSLLAL 98
    FOSB_HUMAN      TSSFVLTCPEVSAFAGAQRTSGSDQPSDPLNSPSLLAL 98
                    ***********************:**************

.. admonition:: Some bioconvert conversions

    :class:`~bioconvert.clustal2fasta.CLUSTAL2FASTA`,
    :class:`~bioconvert.clustal2fasta.CLUSTAL2NEXUS`,
    :class:`~bioconvert.clustal2fasta.CLUSTAL2PHYLIP`,
    :class:`~bioconvert.clustal2fasta.CLUSTAL2STOCKHOLM`,


.. admonition:: Reference

    - https://en.wikipedia.org/wiki/Clustal


.. _format_csv:

CSV
---

:Format: human-readable
:Type: database
:Status: included

A comma-separated values format is a delimited text file that uses a
comma to separate values. See :ref:`format_csv` format page for
details.

.. admonition:: References

    - https://en.wikipedia.org/wiki/Comma-separated_values


.. _format_dsrc:

DSRC
----

:Format: binary
:Status: included
:Type: Compression

DSRC compression dedicated for DNA sequences.

.. admonition:: Bioconvert conversions:

    :class:`~bioconvert.gz2bz2.GZ2BZ2`,
    :class:`~bioconvert.gz2dsrc.GZ2DSRC`
    :class:`~bioconvert.bz22gz.BZ22GZ`,
    :class:`~bioconvert.dsrc2gz.DSRC2GZ`


.. _format_embl:

EMBL
----

:Format: human-readable
:Status: included
:Type: database

EMBL format stores sequence and its annotation together. The start of the
annotation section is marked by a line beginning with the word "ID". The start
of sequence section is marked by a line beginning with the word "SQ".
The "//" (terminator) line also contains no data or comments and designates
the end of an entry.


An example sequence in EMBL format is::

    ID   AB000263 standard; RNA; PRI; 368 BP.
    XX
    AC   AB000263;
    XX
    DE   Homo sapiens mRNA for prepro cortistatin like peptide, complete cds.
    XX
    SQ   Sequence 368 BP;
         acaagatgcc attgtccccc ggcctcctgc tgctgctgct ctccggggcc acggccaccg        60
         ctgccctgcc cctggagggt ggccccaccg gccgagacag cgagcatatg caggaagcgg       120
         caggaataag gaaaagcagc ctcctgactt tcctcgcttg gtggtttgag tggacctccc       180
         aggccagtgc cgggcccctc ataggagagg aagctcggga ggtggccagg cggcaggaag       240
         gcgcaccccc ccagcaatcc gcgcgccggg acagaatgcc ctgcaggaac ttcttctgga       300
         agaccttctc ctcctgcaaa taaaacctca cccatgaatg ctcacgcaag tttaattaca       360
         gacctgaa                                                                368


.. admonition:: Bioconvert conversions:

    :class:`~bioconvert.embl2genbank.EMBL2GENBANK`
    :class:`~bioconvert.genbank2embl.GENBANK2EMBL`

.. admonition:: References

    - ftp://ftp.ebi.ac.uk/pub/databases/embl/release/doc/usrman.txt

.. _format_fam:

FAM
---

:Format: human-readable
:Status: included
:Type: database

The FAM format is used to store sample information accompanying a .bed or
biallelic .pgen binary genotype table. Please see :ref:`format_plink_binary` section.

In brief, it stores the first 6 columns of the PED file. So it is a text file
with no header line, and one line per sample with the following six
fields:

-    Family ID ('FID')
-    Individual ID ('IID'; cannot be '0')
-    Individual ID of father ('0' if father isn't in dataset)
-    Individual ID of mother ('0' if mother isn't in dataset)
-    Sex code ('1' = male, '2' = female, '0' = unknown)
-    Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)

For example::

    1 1000000000 0 0 1 1
    1 1000000001 0 0 1 2

.. _format_faa:

FAA
---

Fasta formatted file storing amino acid sequences. A mutliple protein fasta file
can have the more specific extension mpfa.

.. _format_fasta:

FASTA
-----

:Format: human-readable
:Status: included
:Type: Sequence


FASTA format is one of the most widely used sequence format. It can
stores multiple records of sequence and their identifier.

A sequence entry has a one-line header followed by one or more lines of
sequence. The header must start with the  ">" character. The next word is the
sequence identifier or the accession number; the rest of the line is considered
as description.

The NCBI recommandation do not allowed blank lines in the middle of FASTA files.
Note, however, that some tools can handle blank lines by ignoring them. This is
not recommened to include blank lines though.

There is no standard file extension for a text file containing FASTA formatted sequences. Although
their is a plethora of ad-hoc file extensions: fasta, fas, fa, seq, fsa, fna, ffn, faa, frn, we use only fasta, fa and fst within **Bioconvert** (see :attr:`~bioconvert.core.extensions`). For completeness, *fasta* is the generic fasta file, *fna* stands for fasta nucleic acid, ffn for fasta nucleotide of gene resions, *faa* for fasta amino acid, *frn* for fasta non-coding RNA, etc.

An example sequence in FASTA format is::

    >X65923.1 H.sapiens fau mRNA
    TTCCTCTTTCTCGACTCCATCTTCGCGGTAGCTGGGACCGCCGTTCAGTCGCCAATATGCAGCTCTTTGT
    CCGCGCCCAGGAGCTACACACCTTCGAGGTGACCGGCCAGGAAACGGTCGCCCAGATCAAGGCTCATGTA
    GCCTCACTGGAGGGCATTGCCCCGGAAGATCAAGTCGTGCTCCTGGCAGGCGCGCCCCTGGAGGATGAGG
    CCACTCTGGGCCAGTGCGGGGTGGAGGCCCTGACTACCCTGGAAGTAGCAGGCCGCATGCTTGGAGGTAA
    AGTTCATGGTTCCCTGGCCCGTGCTGGAAAAGTGAGAGGTCAGACTCCTAAGGTGGCCAAACAGGAGAAG
    AAGAAGAAGAAGACAGGTCGGGCTAAGCGGCGGATGCAGTACAACCGGCGCTTTGTCAACGTTGTGCCCA
    CCTTTGGCAAGAAGAAGGGCCCCAATGCCAACTCTTAAGTCTTTTGTAATTCTGGCTTTCTCTAATAAAA
    AAGCCACTTAGTTCAGTCAAAAAAAAAA

In this example, the header (also known as description line) is formatted as::

    >ID description

Many variants of FASTA formats exists but differ only in the way the header is
written. All starts with the ">" sign though. We can cite a few variants here
below (for simplicity we give only puit 2 lines per sequence).


The **NCBI style** defines the identifier with database name, entry ID and
optional accession or sequence version number separated by pipes::

    >embl|X65923|X65923.1 H.sampiens fau mRNA
    TTCCTCTTTCTCGACTCCATCTTCGCGGTAGCTGGGACCGCCGTTCAGTCGCCAATATGCAGCTCTTTGT
    CCGCGCCCAGGAGCTACACACCTTCGAGGTGACCGGCCAGGAAACGGTCGCCCAGATCAAGGCTCATGTA

List of NCBI FASTA database are listed in https://tinyurl.com/y6wrzyad

The **GI style** is the same as NCBI style except that the sequence GI code is
given instead of the entry ID::

    >gi|31302|gnl|genbank|X65923 (X65923.1) H.sampiens fau mRNA
    TTCCTCTTTCTCGACTCCATCTTCGCGGTAGCTGGGACCGCCGTTCAGTCGCCAATATGCAGCTCTTTGT
    CCGCGCCCAGGAGCTACACACCTTCGAGGTGACCGGCCAGGAAACGGTCGCCCAGATCAAGGCTCATGTA

There is also a **CGC-style** FASTA format (not to be confused with
the :ref:`format_gcg` format). Its header includes an optional database
name as part of the identifier by  using the : sign::

      >DATABASE_NAME:DI accession description

      >embl:X65923 X65923.1 H.sapiens fau mRNA
      TTCCTCTTTCTCGACTCCATCTTCGCGGTAGCTGGGACCGCCGTTCAGTCGCCAATATGCAGCTCTTTGT
      CCGCGCCCAGGAGCTACACACCTTCGAGGTGACCGGCCAGGAAACGGTCGCCCAGATCAAGGCTCATGTA

And more generally, we have the FASTA with accession and description style.
The accession number or sequence version included after the identifier::

    >X65923 X65923.1 H.sapiens fau mRNA
    TTCCTCTTTCTCGACTCCATCTTCGCGGTAGCTGGGACCGCCGTTCAGTCGCCAATATGCAGCTCTTTGT
    CCGCGCCCAGGAGCTACACACCTTCGAGGTGACCGGCCAGGAAACGGTCGCCCAGATCAAGGCTCATGTA
    GCCTCACTGGAGGGCATTGCCCCGGAAGATCAAGTCGTGCTCCTGGCAGGCGCGCCCCTGGAGGATGAGG
    CCACTCTGGGCCAGTGCGGGGTGGAGGCCCTGACTACCCTGGAAGTAGCAGGCCGCATGCTTGGAGGTAA
    AGTTCATGGTTCCCTGGCCCGTGCTGGAAAAGTGAGAGGTCAGACTCCTAAGGTGGCCAAACAGGAGAAG
    AAGAAGAAGAAGACAGGTCGGGCTAAGCGGCGGATGCAGTACAACCGGCGCTTTGTCAACGTTGTGCCCA
    CCTTTGGCAAGAAGAAGGGCCCCAATGCCAACTCTTAAGTCTTTTGTAATTCTGGCTTTCTCTAATAAAA
    AAGCCACTTAGTTCAGTCAAAAAAAAAA

.. note:: original FASTA format may include comments with the ; sign. This
   is not supported anymore in most programs.

.. admonition:: Bioconvert conversions

    :class:`~bioconvert.fastq2fasta.FASTQ2FASTA`, :class:`~bioconvert.fasta2fasta.FASTA2FASTQ`,
    :class:`~bioconvert.fasta2clustal.FASTA2CLUSTAL`, :class:`~bioconvert.fasta2nexus.FASTA2NEXUS`,
    :class:`~bioconvert.fasta2twobit.FASTA2TWOBIT`

.. seealso:: :ref:`format_fastq` and :ref:`format_qual`
.. admonition::  References

    -  http://en.wikipedia.org/wiki/FASTA_format
    -  NCBI recommandations: https://tinyurl.com/y6wrzyad


.. _format_fastg:

FastG
-----

:Format:
:Status: not included
:Type: assembly


FastG is a Graph format used to faithfully representing genome
assemblies in the face of allelic polymorphism and assembly uncertainty.

:reference: http://fastg.sourceforge.net/FASTG_Spec_v1.00.pdf


.. _format_fastq:

FastQ
-----

:Format: human-readable
:Status: included
:Type: Sequence


FASTQ is a text-based format for storing both biological sequence (usually
nucleotide sequence) and its corresponding quality scores (:ref:`format_qual`).
A FASTQ format can contain several sequences. All FASTQ variations are in the
formatting of the quality scores. Currently, the recommended variant is the
Sanger encoding also used by Illumina 1.8. It encodes the Phred quality score
from 0 to 93 using ASCII 33 to 126. This format is also refered
to PHRED+33 meaning there is an offset of 33 in the ASCII code. Other variants
such as FASTQ-solexa or earlier Illumina versions. Currently conversions included
in **Bioconvert** do not need to be aware of the quality score encoding.


A FASTQ file uses four lines per sequence:

1. a '@' character, followed by a sequence identifier and an optional description
2. the raw sequence letters.
3. a '+' character, optionally followed by the same sequence identifier (and any description)
4. quality values for the sequence in Line 2

An example sequence in FASTQ format is::

    @SEQUENCE_ID1
    GTGGAAGTTCTTAGGGCATGGCAAAGAGT
    +
    FAFFADEDGDBGEGGBCGGHE>EEBA@@=
    @SEQUENCE_ID2
    GTGGAAGTTCTTAGG
    +
    FAFFADEDGDBGEGG


.. admonition:: Bioconvert conversions

    :class:`~bioconvert.fastq2fasta.FASTQ2FASTA`, :class:`~bioconvert.fasta2fasta.FASTA2FASTQ`

.. seealso:: :ref:`format_fasta` and :ref:`format_qual`


.. admonition:: References

    - https://en.wikipedia.org/wiki/FASTQ_format

.. _format_genbank:

GENBANK
-------

:Format: human-readable
:Status: included
:Type: annotation/sequence

GenBank format (GenBank Flat File Format) stores sequence and its annotation
together. The start of the annotation section is marked by a line beginning with
the word *LOCUS*. The start of sequence section is marked by a line beginning
with the word *ORIGIN* and the end of the section is marked by a line with only
"//".

GenBank format for protein has been renamed **GenPept**.

An example sequence in GenBank format is::

    LOCUS       AB000263                 368 bp    mRNA    linear   PRI 05-FEB-1999
    DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete
                cds.
    ACCESSION   AB000263
    ORIGIN
            1 acaagatgcc attgtccccc ggcctcctgc tgctgctgct ctccggggcc acggccaccg
           61 ctgccctgcc cctggagggt ggccccaccg gccgagacag cgagcatatg caggaagcgg
          121 caggaataag gaaaagcagc ctcctgactt tcctcgcttg gtggtttgag tggacctccc
          181 aggccagtgc cgggcccctc ataggagagg aagctcggga ggtggccagg cggcaggaag
          241 gcgcaccccc ccagcaatcc gcgcgccggg acagaatgcc ctgcaggaac ttcttctgga
          301 agaccttctc ctcctgcaaa taaaacctca cccatgaatg ctcacgcaag tttaattaca
          361 gacctgaa
    //



.. admonition:: Bioconvert conversions

    :class:`~bioconvert.genbank2fasta.GENBANK2FASTA`,
    :class:`~bioconvert.genbank2embl.GENBANK2EMBL`

.. admonition:: References:

    - https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html


.. _format_genpept:

GENPEPT
-------

see :ref:`format_genbank`


.. _gfa_format:

GFA
---

:Format: human-readable
:Status: included
:Type: assembly graph

The Graphical Fragment Assembly (GFA) can be used to represent genome
assemblies. GFA stores sequence graphs as the product of an
assembly, a representation of variation in genomes, splice graphs in genes, or
even overlap between reads from long-read sequencing technology.

The GFA format is a tab-delimited text format for describing a set of sequences
and their overlap. The first field of the line identifies the type of the line.
**Header** lines start with H. **Segment** lines start with S. **Link** lines start with L.
A **containment** line starts with C. A **path** line starts with P.

- Segment a continuous sequence or subsequence.
- Link an overlap between two segments. Each link is from the end of one segment to the beginning of another segment. The link stores the orientation of each segment and the amount of basepairs overlapping.
- Containment an overlap between two segments where one is contained in the other.
- Path an ordered list of oriented segments, where each consecutive pair of oriented segments are supported by a link record.

See details in the reference above.

Example:

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


GFA2 is a generalization of GFA that allows one to specify an assembly graph in
either less detail, e.g. just the topology of the graph, or more detail, e.g.
the multi-alignment of reads giving rise to each sequence. It is further
designed to be a able to represent a string graph at any stage of assembly, from
the graph of all overlaps, to a final resolved assembly of contig paths with
multi-alignments. Apart from meeting these needs, the extensions also supports
other assembly and variation graph types.

Like GFA, GFA2 is tab-delimited in that every lexical token is separated from
the next by a single tab.

.. admonition:: Bioconvert conversions

    :class:`~bioconvert.gfa2fasta.GFA2FASTA`

.. admonition:: References:

    - http://gfa-spec.github.io/GFA-spec/,

.. _format_gff:

GTF
---

:Format: human-readable
:Status: included
:Type: Annotation


GTF2 (General Feature Format version 2) is a file format used to represent genomic features and their locations in a genome. It is a tab-delimited text file that contains one line for each genomic feature, with each line consisting of nine fields separated by tabs.

The fields in a GTF2 file are as follows:

* Seqid: The identifier of the genomic sequence.
* Source: The source of the annotation.
* Feature: The type of feature.
* Start: The starting position of the feature.
* End: The ending position of the feature.
* Score: A score associated with the feature.
* Strand: The strand on which the feature is located.
* Phase: The phase of the feature, if applicable.
* Attributes: A list of attributes associated with the feature, encoded as a semicolon-separated list of key-value pairs.

GFF
---

:Format: human-readable
:Status: included
:Type: Annotation

GFF is a standard file format for storing genomic features in a text file. GFF
stands for Generic Feature Format. It is 9 column tab-delimited file, each line of which corresponds to an annotation, or feature.

The GFF v2 is deprecated and v3 should be used instead. In particular, GFF2 is sunable to deal with the three-level hierarchy of gene -> transcript -> exon.

The first line is a comment (starting with #) followed by a series of data lines, each of which correspond to an annotation. Here is an example::

    ##gff-version 3
    ctg123  .  exon  1300  1500  .  +  .  ID=exon00001
    ctg123  .  exon  1050  1500  .  +  .  ID=exon00002
    ctg123  .  exon  3000  3902  .  +  .  ID=exon00003
    ctg123  .  exon  5000  5500  .  +  .  ID=exon00004
    ctg123  .  exon  7000  9000  .  +  .  ID=exon00005

The header is compulsary and following lines must have 9 columns as follows:

1. **seqname** - The name of the sequence (e.g. chromosome) on which the feature
   exists. Any string can be used. For example, *chr1*, *III*, *contig1112.23*.
   Any character not in  ``[a-zA-Z0-9.:^*$@!+_?-|]`` must be escaped with the %
   character followed by its hexadecimal value.
2. **source** - The source of this feature. This field will normally be used
   to indicate the program making the prediction, or if it comes from public
   database annotation, or is experimentally verified, etc. If there is no
   source, use the . character.
3. **feature** - The feature type name. Equivalent to BED’s name field.  For example, *exon*, etc. Should be a term from the lite sequence ontology (SOFA).
4. **start** - The one-based starting position of feature on seqname.
   bedtools uses a one-based position and BED uses a zero-based start position.
5. **end** - The one-based ending position of feature on seqname.
6. **score** - A score assigned to the GFF feature.
7. **strand** - Defines the strand. Use +, - or .
8. **frame/phase** - The frame of the coding sequence. Use 0, 1, 2. The phase
   is one
   of the integers 0, 1, or 2, indicating the number of bases that should be
   removed from the beginning of this feature to reach the first base of the
   next codon.
9. **attribute** - A list of feature attributes in the format tag=value
   separated by semi columns.
   All non-printing characters in such free text value strings (e.g. newlines,
   tabs, control characters, etc) must be explicitly represented by their
   C (UNIX) style backslash-escaped representation (e.g. newlines as ‘n’,
   tabs as ‘t’). Tabs must be replaced with %09 URL escape. There are predefined
   tags:

   - ID: unique identifier of the feature.
   - Name: name of the feature
   - Alias
   - Parent: can be used to group exons into transcripts, transcripts into
     genes and so on.
   - Target
   - Gap
   - Derives_from
   - Note
   - Dbxref
   - Ontology_term

   Multiple attributes of the same type are separated by comma.
   Case sensitive: Parent is difference from parent.

.. admonition:: Bioconvert conversions:

    - :class:`~bioconvert.gff22gff3.GFF22GFF3`,
      :class:`~bioconvert.gff32gff2.GFF32GFF2`


.. admonition:: References:

    - http://gmod.org/wiki/GFF2
    - http://gmod.org/wiki/GFF3
    - http://www.sanger.ac.uk/resources/software/gff/spec.html
    - https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md



.. _format_gz:

GZ
--

:Format: binary
:Status: included
:Type: Compression

**gzip** is a file compression program that is based on the DEFLATE algorithm, which is a combination of LZ77 and Hufmfman coding.

.. admonition:: Bioconvert conversions:

    :class:`~bioconvert.gz2bz2.GZ2BZ2`,
    :class:`~bioconvert.gz2dsrc.GZ2DSRC`
    :class:`~bioconvert.bz22gz.BZ22GZ`,
    :class:`~bioconvert.dsrc2gz.DSRC2GZ`


.. _format_json:

JSON
----

:Format: human-readable
:Status: included
:Type: database

JSON format stands for Javascript Object Notation. Basic data types used in JSON:

- **Number**: a signed decimal number that may contain a fractional part and may use
  exponential E notation, but cannot include non-numbers such as NaN. The format
  makes no distinction between integer and floating-point. JavaScript uses a
  double-precision floating-point format for all its numeric values, but other
  languages implementing JSON may encode numbers differently.
- **String**: a sequence of zero or more Unicode characters. Strings are delimited
  with double-quotation marks and support a backslash escaping syntax.
- **Boolean**: either of the values true or false
- **Array**: an ordered list of zero or more values, each of which may be of any type.
  Arrays use square bracket notation and elements are comma-separated.
- **Object**: an unordered collection of name–value pairs where the names (also called
  keys) are strings. Since objects are intended to represent associative
  arrays, it is recommended that each key is unique
  within an object. Objects are delimited with curly brackets and use commas to
  separate each pair, while within each pair the colon ':' character separates the
  key or name from its value.
- **null**: An empty value, using the word null

Limited whitespace is allowed and ignored around or between syntactic elements
(values and punctuation, but not within a string value). Only four specific
characters are considered whitespace for this purpose: space, horizontal tab,
line feed, and carriage return. In particular, the byte order mark must not be
generated by a conforming implementation (though it may be accepted when parsing
JSON). JSON does not provide syntax for comments.

Example::

    {
    "database": "AB",
    "date": "13-10-2010",
    "entries":
        [
          {
            "ID": 1,
            "coverage": 10
          },
          {
            "ID": 2,
            "coverage": 15
          }
        ]
    }



.. admonition:: Bioconvert conversions

    :class:`~bioconvert.json2yaml.JSON2YAML`,
    :class:`~bioconvert.yaml2json.YAML2JSON`.


.. admonition:: References

    - https://en.wikipedia.org/wiki/JSON

.. _format_maf_mutation:

MAF (Mutation Annotation Format)
--------------------------------

:Format: human-readable
:Status: not included
:Type: multiple alignement


.. admonition:: reference:

    - https://software.broadinstitute.org/software/igv/MutationAnnotationFormat


.. _format_maf:

MAF (Multiple Alignement Format)
--------------------------------

:Format: human-readable
:Status: included
:Type: phylogeny

The Multiple Alignment Format stores a series of multiple alignments.

.. warning:: Not to be confused with :ref:`format_maf_mutation`

Here are some rules about the MAF syntax:

- It is line-oriented.
- Each multiple alignment ends with a blank line.
- Each sequence in an alignment is on a single line, which can get quite
  long, but there is no length limit.
- Words in a line are delimited by any white space.
- Lines starting with # are considered to be comments.
- Lines starting with ## can be ignored by most programs, but contain meta-data of one
  form or another.
- The file is divided into paragraphs that terminate in a blank line.
- Within a paragraph, the first word of a line indicates its type.

Each multiple alignment is in a separate paragraph that begins with an **a** line and contains an **s** line for each sequence in the multiple alignment.

Some MAF files may contain other optional line types:

- **i** line contains information about what is in the aligned
  species DNA before and after the immediately preceding **s** line
- **e** line contains information about the size of the gap
  between the alignments that span the current block
- **q** line indicates the quality of each aligned base for the species.

Here is an example of **s** lines (alignment block)::

    s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca
    s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
    s baboon         249182 13 +   4622798 gcagctgaaaaca
    s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA

The **s** and **a** lines define a multiple alignment. The columns of the **s** lines
have the following fields:

- **src**:  The name of one of the source sequences for the alignment. The form 'database.chromosome' allows automatic creation of links to other assemblies in some browsers.
- **start**: The start of the aligning region in the source sequence. This is a zero-based number. If the strand field is "-" then this is the start relative to the reverse-complemented source sequence (see Coordinate Transforms).
- **size**:  The size of the aligning region in the source sequence. This number is equal to the number of non-dash characters in the alignment text field below.
- **strand**: Either + or -. If -, then the alignment is to the reverse-complemented source.
- **srcSize**: The size of the entire source sequence, not just the parts involved in the alignment.
- **text**: The nucleotides (or amino acids) in the alignment and any insertions (dashes).

Lines starting with **i** give information about what's happening before and after this block in the aligning species::

    s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca
    s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
    i panTro1.chr6 N 0 C 0
    s baboon         249182 13 +   4622798 gcagctgaaaaca
    i baboon       I 234 n 19

The **i** lines contain information about the context of the sequence lines immediately preceding them. The following fields are defined by position rather than name=value pairs:

- **src**: The name of the source sequence for the alignment. Should be the
  same as the **s** line immediately above this line.
- **leftStatus**: A character that specifies the relationship between the sequence
  in this block and the sequence that appears in the previous block.
- **leftCount**: Usually the number of bases in the aligning species between the
  start of this alignment and the end of the previous one.
- **rightStatus**: A character that specifies the relationship between the sequence
  in this block and the sequence that appears in the subsequent block.
- **rightCount**: Usually the number of bases in the aligning species between the
  end of this alignment and the start of the next one.

The status characters can be one of the following values::

    C: the sequence before or after is contiguous with this block.
    I: there are bases between the bases in this block and the one before or
       after it.
    N: this is the first sequence from this src chrom or scaffold.
    n: this is the first sequence from this src chrom or scaffold but it is
       bridged by another alignment from a different chrom or scaffold.
    M: there is missing data before or after this block (Ns in the sequence).
    T: the sequence in this block has been used before in a previous block
       (likely a tandem duplication)

Lines starting with **e** gives information about empty parts of the alignment
block::

    s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca
    e mm4.chr6     53310102 13 + 151104725 I

The **e** lines indicate that there isn't aligning DNA for a species but that
the current block is bridged by a chain that connects blocks before and after
this block. The following fields are defined by position rather than name=value pairs.

- **src**: The name of one of the source sequences for the alignment.
- **start**: The start of the non-aligning region in the source sequence.
  This is a zero-based number. If the strand field is "-" then this
  is the start relative to the reverse-complemented source sequence
  (see Coordinate Transforms).
- **size**: The size in base pairs of the non-aligning region in the
  source sequence.
- **strand**: Either + or -. If -, then the alignment is to the reverse-complemented source.
- **srcSize**: The size of the entire source sequence, not just the parts involved
  in the alignment; alignment and any insertions (dashes) as well.
- *status**: A character that specifies the relationship between the non-aligning
  sequence in this block and the sequence that appears in the previous
  and subsequent blocks.

The status character can be one of the following values::

    C: the sequence before and after is contiguous implying that this region
       was either deleted in the source or inserted in the reference sequence.
       The browser draws a single line or a "-" in base mode in these blocks.
    I: there are non-aligning bases in the source species between chained alignment
       blocks before and after this block. The browser shows a double line
       or "=" in base mode.
    M: there are non-aligning bases in the source and more than 90% of them are Ns in
       the source. The browser shows a pale yellow bar.
    n: there are non-aligning bases in the source and the next aligning block starts
       in a new chromosome or scaffold that is bridged by a chain between still
       other blocks. The browser shows either a single line or a double line based
       on how many bases are in the gap between the bridging alignments.

Lines starting with **q** -- information about the quality of each aligned base for the species::

    s hg18.chr1                  32741 26 + 247249719 TTTTTGAAAAACAAACAACAAGTTGG
    s panTro2.chrUn            9697231 26 +  58616431 TTTTTGAAAAACAAACAACAAGTTGG
    q panTro2.chrUn                                   99999999999999999999999999
    s dasNov1.scaffold_179265     1474  7 +      4584 TT----------AAGCA---------
    q dasNov1.scaffold_179265                         99----------32239---------

The **q** lines contain a compressed version of the actual raw quality data, representing
the quality of each aligned base for the species with a single character of 0-9 or F.
The following fields are defined by position rather than name=value pairs:

- **src**: The name of the source sequence for the alignment. Should be the same as the "s" line immediately preceding this line.
- **value**: A MAF quality value corresponding to the aligning nucleotide acid in
  the preceding "s" line. Insertions (dashes) in the preceding "s" line are represented
  by dashes in the "q" line as well. The quality value can be "F" (finished sequence)
  or a number derived from the actual quality scores (which range from 0-97) or the
  manually assigned score of 98. These numeric values are calculated as::

    MAF quality value = min( floor(actual quality value/5), 9 )

This results in the following mapping::

    MAF quality value     Raw quality score range     Quality level
    0-8     0-44     Low
    9     45-97     High
    0     98     Manually assigned
    F     99     Finished

A Simple Example (three alignment blocks derived from five starting sequences).
Repeats are shown as lowercase, and each block may have a subset of the input sequences.
All sequence columns and rows must contain at least one nucleotide (no columns or rows that contain only insertions)::

    ##maf version=1 scoring=tba.v8
    # tba.v8 (((human chimp) baboon) (mouse rat))

    a score=23262.0
    s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
    s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
    s baboon         116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
    s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
    s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG

    a score=5062.0
    s hg18.chr7    27699739 6 + 158545518 TAAAGA
    s panTro1.chr6 28862317 6 + 161576975 TAAAGA
    s baboon         241163 6 +   4622798 TAAAGA
    s mm4.chr6     53303881 6 + 151104725 TAAAGA
    s rn3.chr4     81444246 6 + 187371129 taagga

    a score=6636.0
    s hg18.chr7    27707221 13 + 158545518 gcagctgaaaaca
    s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
    s baboon         249182 13 +   4622798 gcagctgaaaaca
    s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA


.. admonition:: References

    - https://github.com/peterjc/maf2sam/
    - https://github.com/arq5x/nanopore-scripts/master/maf-convert.py
    - Example and doc from https://genome.ucsc.edu/FAQ/FAQformat.html#format5


.. _format_map:

MAP
---

:Format: human-readable
:Status: included
:Type: Genotypic

PLINK is a very widely used application for analyzing genotypic data.

The fields in a MAP file are:

- Chromosome
- Marker ID
- Genetic distance
- Physical position

Example of a MAP file of the standard PLINK format::

    21     rs11511647   0          26765
    X      rs3883674    0           32380
    X      rs12218882   0           48172
    9      rs10904045   0           48426
    9      rs10751931   0           49949
    8      rs11252127   0           52087
    10     rs12775203   0           52277
    8      rs12255619   0           52481

.. admonition:: References

    - http://www.gwaspi.org/?page_id=145

.. _format_newick:

NEWICK
------

:Format: human-readable
:Status: included
:Type: phylogeny

Newick format is typically used for tools like PHYLIP and is a minimal
definition for a phylogenetic tree. It is a way of representing
graph-theoretical trees with edge lengths using parentheses and commas.


.. image:: _static/NewickExample.svg

::


    (,,(,));                              no nodes are named
    (A,B,(C,D));                          leaf nodes are named
    (A,B,(C,D)E)F;                        all nodes are named
    (:0.1,:0.2,(:0.3,:0.4):0.5);          all but root node have a distance to parent
    (:0.1,:0.2,(:0.3,:0.4):0.5):0.0;      all have a distance to parent
    (A:0.1,B:0.2,(C:0.3,D:0.4):0.5);      distances and leaf names (popular)
    (A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;    distances and all names
    ((B:0.2,(C:0.3,D:0.4)E:0.5)A:0.1)F;   a tree rooted on a leaf node (rare)

.. admonition:: Bioconvert conversions

    :class:`~bioconvert.newick2nexus.NEWICK2NEXUS`,
    :class:`~bioconvert.newick2phyloxml.NEWICK2PHYLOXML`


.. admonition:: References

    - https://en.wikipedia.org/wiki/Newick_format
    - http://evolution.genetics.washington.edu/phylip/newicktree.html

.. _format_nexus:

NEXUS
-----

:Format: human-readable
:Status: included
:Type: phylogeny

The NEXUS multiple alignment format, also known as PAUP format is used to
multiple alignment or phylogentic trees.


After a header to indicate the format (#NEXUS ), blocks are stored
and start with *Begin NAME;* and end with *END;*

Example of a DNA alignment::

    #NEXUS
    Begin data;
    Dimensions ntax=4 nchar=15;
    Format datatype=dna missing=? gap=-;
    Matrix
    Species1   atgctagctagctcg
    Species2   atgcta??tag-tag
    Species3   atgttagctag-tgg
    Species4   atgttagctag-tag
    ;
    End;

It can be used to store phylogenetic trees using the TREES block::

    #NEXUS
    BEGIN TAXA;
      TAXLABELS A B C;
    END;

    BEGIN TREES;
      TREE tree1 = ((A,B),C);
    END;



.. admonition:: Bioconvert conversions

    :class:`~bioconvert.nexus2clustal.NEXUS2CLUSTAL`,
    :class:`~bioconvert.nexus2newick.NEXUS2NEWICK`,
    :class:`~bioconvert.nexus2phylip.NEXUS2PHYLIP`,
    :class:`~bioconvert.nexus2phyloxml.NEXUS2PHYLIPXML`,


.. admonition:: References

    - https://en.wikipedia.org/wiki/Nexus_file



.. _format_ods:

ODS
---

:Format: human-readable
:Status: included
:Type: Sequence

ODS stands for OpenDocument Spreadsheet (.ods) file format. It should be
equivalent to the :ref:`format_xls` format.


.. admonition:: Bioconvert conversions

    :class:`~bioconvert.json2yaml.JSON2YAML`,
    :class:`~bioconvert.yaml2json.YAML2JSON`.


.. _format_paf:

PAF (Pairwise mApping Format)
-----------------------------

:Format: human-readable
:Status: included
:Type: mapping


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

.. admonition:: Bioconvert conversion

    :class:`~bioconvert.sam2paf.SAM2PAF`


.. admonition:: References:

    - https://github.com/lh3/miniasm/blob/master/PAF.md

.. _format_pdb:

PDB
---

.. todo:: coming soon


.. _format_ped:

PED
---

:Format: human-readable
:Status: included
:Type: Genotypic

PLINK is a very widely used application for analyzing genotypic data.

The fields in a PED file are:

- Family ID
- Sample ID
- Paternal ID
- Maternal ID
- Sex (1=male; 2=female; other=unknown)
- Affection (0=unknown; 1=unaffected; 2=affected)
- Genotypes (space or tab separated, 2 for each marker. 0=missing)


Example of a PED file of the standard PLINK format::

    FAM1    NA06985 0   0   1   1   A   T   T   T   G   G   C   C   A   T   T   T   G   G   C   C
    FAM1    NA06991 0   0   1   1   C   T   T   T   G   G   C   C   C   T   T   T   G   G   C   C
    0       NA06993 0   0   1   1   C   T   T   T   G   G   C   T   C   T   T   T   G   G   C   T
    0       NA06994 0   0   1   1   C   T   T   T   G   G   C   C   C   T   T   T   G   G   C   C
    0       NA07000 0   0   2   1   C   T   T   T   G   G   C   T   C   T   T   T   G   G   C   T
    0       NA07019 0   0   1   1   C   T   T   T   G   G   C   C   C   T   T   T   G   G   C   C
    0       NA07022 0   0   2   1   C   T   T   T   G   G   0   0   C   T   T   T   G   G   0   0
    0       NA07029 0   0   1   1   C   T   T   T   G   G   C   C   C   T   T   T   G   G   C   C
    FAM2    NA07056 0   0   0   2   C   T   T   T   A   G   C   T   C   T   T   T   A   G   C   T
    FAM2    NA07345 0   0   1   1   C   T   T   T   G   G   C   C   C   T   T   T   G   G   C   C



.. _format_phyloxml:

PHYLOXML
--------

:Format: human-readable
:Status: included
:Type: phylogeny


PhyloXML is an XML language for the analysis, exchange, and storage of phylogenetic
trees.

A shortcoming of formats such as Nexus and Newick is a lack of a
standardized means to annotate tree nodes and branches with distinct data fields
(species names, branch lengths, multiple support values). A well defined XML format
addresses these problems in a general and extensible manner and allows for
interoperability between specialized and general purpose software.

Here is an example (source https://en.wikipedia.org/wiki/PhyloXML) ::


    <phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd"
    xmlns="http://www.phyloxml.org">
    <phylogeny rooted="true">
      <name>example from Prof. Joe Felsenstein's book "Inferring Phylogenies"</name>
      <description>MrBayes based on MAFFT alignment</description>
      <clade>
         <clade branch_length="0.06">
            <confidence type="probability">0.88</confidence>
            <clade branch_length="0.102">
               <name>A</name>
            </clade>
            <clade branch_length="0.23">
               <name>B</name>
            </clade>
          </clade>
          <clade branch_length="0.5">
            <name>C</name>
          </clade>
        </clade>
      </phylogeny>
    </phyloxml>


.. admonition:: Bioconvert conversions

    :class:`~bioconvert.phyloxml2nexus.PHYLOXML2NEXUS`
    :class:`~bioconvert.phyloxml2newick.PHYLOXML2NEWICK`


.. admonition:: References

    - http://www.phyloxml.org/
    - https://en.wikipedia.org/wiki/PhyloXML


.. _format_phylip:

PHYLIP
------

:Format: human-readable
:Status: included
:Type: phylogeny / alignement

The PHYLIP format stores a multiple sequence alignement.

It is a plain test format with a header describing the dimensions of the
alignment followed by the mutliple sequence alignment. The following sequence
is exactly 10 characters long (padded wit spaces if needed).

PHYLIP does not support blank lines between header and the alignment.

In the header, the first integer defines the number of sequences.
The second intefer defines the number of alignments. There are several
spaces between the two integers.

Here is an example::

       5   50
    Seq0000  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC
    Seq0001  ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
    Seq0002  TGTCGGACCT AAGTATTGAG TACAACGGTG TATTCCAGCG GTGGAGAGGT
    Seq0003  CTATTTTTCC GGTTGAAGGA CTCTAGAGCT GTAAAGGGTA TGGCCATGTG
    Seq0004  CTAAGCGCGG GCGGATTGCT GTTGGAGCAA GGTTAAATAC TCGGCAATGC


.. admonition:: Bioconvert conversions

    :class:`~bioconvert.phylip2clustal.PHYLIP2CLUSTAL`,
    :class:`~bioconvert.phylip2fasta.PHYLIP2FASTA`,
    :class:`~bioconvert.phylip2nexus.PHYLIP2NEXUS`,
    :class:`~bioconvert.phylip2stockholm.PHYLIP2STOCKHOLM`


.. admonition:: References

    - http://www.phyloxml.org/
    - https://en.wikipedia.org/wiki/PhyloXML

.. _format_plink_flat:

PLINK flat files (MAP/PED)
--------------------------

:Format: human-readable
:Status: included
:Type: genotypic

PLINK is a used application for analyzing genotypic data. It can be considered
the de-facto standard of the field.

The standard PLINK files can be a bundle of plain text files (PED & MAP dataset,
or its transpose, TPED & :ref:`format_fam` dataset), or a bundle of binary
files (BED, :ref:`format_bim` & :ref:`format_fam`) as explained in :ref:`format_plink_binary`.

PLINK provides commands to convert between text and binary formats. In
Bioconvert, you can use the **plink2bpblink** conversion::

    bioconvert plink2bplink input_prefix output_prefix

.. note:: Since there are several input and output files, we do not provide the
   extension. Instead, we use the prefix filename.


Since PLINK files do not specify for a variant which allele is reference and which is
alternative, importing data to a variant tools project requires matching each
variant to the reference sequence to determine reference and alternative
alleles


The Genotypic data are separated in two flat files: MAP and PED.

The MAP files describes the SNPs and contains those fields:

- chromosome number (integer)
- SNP marker ID (string)
- SNP generit position (cM) (float)
- SNP physical position (bp)

It is spaced or tabulated file with 4 columns. All SNPs must be ordered by physical
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

Then, additional columns can be:

- columns 7 and 8 code for the observed alleles at SNP1
- comumns 9 and 10 code for the observed alleles at SNP2 and so on

Missing data are coded as "0 0". So we have N lines times 2L + 6 columns where N is
the number of individuals and L the numbers of SNPs


Given a .ped file (plink format), we can convert it into the 012 format (0
hom ancestral), 1 het, 2 dom derived using ::

    plink --file [.ped/.map fileset prefix] --recodeA --out [output prefix]

.. admonition:: Bioconvert conversions

    :class:`~bioconvert.plink2bplink.PLINK2BPLINK`,
    :class:`~bioconvert.bplink2plink.BPLINK2PLINK`,
    :class:`~bioconvert.plink2vcf.PLINK2VCF`

.. admonition:: References

    - http://www.gwaspi.org/?page_id=145
    - https://vatlab.github.io/


.. _format_plink_binary:

PLINK binary files (BED/BIM/FAM)
--------------------------------

:Format: human-readable and biarny
:Status: included
:Type: genotypic

PLINK binary format (BED, :ref:`format_bim` and :ref:`format_fam`) is a valid input for many software. If you have the :ref:`format_plink_flat` version, use PLINK to convert text to binary format if necessary. In Bioconvert, you can use the **plink2bpblink** as explained in the :ref:`format_plink_flat` section.

Here, the BED file is binary and is not to be confused with the BEDGRAPH format.

.. admonition:: Bioconvert conversions

    :class:`~bioconvert.plink2bplink.PLINK2BPLINK`,
    :class:`~bioconvert.bplink2plink.BPLINK2PLINK`,
    :class:`~bioconvert.plink2vcf.PLINK2VCF`


.. admonition:: Reference

    - http://zzz.bwh.harvard.edu/plink/tutorial.shtml

.. _format_qual:

QUAL
----

:Format: human-readable
:Status: included
:Type: Sequence

QUAL files include qualities of each nucleotide in :ref:`format_fasta` format.

.. admonition:: Bioconvert conversions

    - :class:`~bioconvert.fastq2fasta.FASTQ2FASTA`
    - :class:`~bioconvert.fasta2fasta.FASTA2FASTQ`

.. seealso:: :ref:`format_fasta` and :ref:`format_fastq`



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


.. admonition:: Bioconvert conversions

    :class:`~bioconvert.bam2sam.BAM2SAM`, :class:`~bioconvert.sam2bam.SAM2BAM`

.. admonition::  References

    - http://samtools.github.io/hts-specs/SAMv1.pdf
    - http://genome.ucsc.edu/goldenPath/help/bam.html


.. _format_scf:

SCF
---

:Format: human readable
:Status: included
:Type: alignment


Trace File Format - Sequence Chromatogram Format (SCF) is a binary file
containing raw data output from automated sequencing instruments.

This converter was translated from BioPerl.


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

.. admonition:: Bioconvert conversions

    :class:`~bioconvert.convert.scf2fastq.SCF2FASTQ`,
    :class:`~bioconvert.convert.scf2fasta.SCF2FASTA`.

.. admonition:: References

    - https://wiki.nci.nih.gov/display/TCGA/Sequence+trace+files
    - http://staden.sourceforge.net/manual/formats_unix_2.html


.. _format_sra:

SRA
---

The Sequence Read Archive (SRA) makes biological sequence data available to the
research community. It stores raw sequencing data and alignment
information from high-throughput sequencing platforms, including Roche 454 GS
System, Illumina Genome Analyzer, Applied Biosystems SOLiD System, Helicos
Heliscope, Complete Genomics, and Pacific Biosciences SMRT.

It is not a format per se but is included in Bioconvert by allowing the
retrieval of sequencing data given a SRA identifier::

    bioconvert sra2fastq <SRA_ID>

This will retrieve the fastq reads (single read or paired end data).


.. admonition:: Bioconvert conversions

    :class:`~bioconvert.sra2fastq.SRA2FASTQ`

.. admonition:: Reference:

    - https://www.ncbi.nlm.nih.gov/sra

.. _format_tsv:

TSV
---

:Format: human readable
:Type: database
:Status: included

A tab-separated values format is a delimited text file that uses a
tab character to separate values. See :ref:`format_csv` format page for
details.


.. admonition:: Bioconvert conversions:

    :class:`~bioconvert.tsv2csv.TSV2CSV`,



.. admonition:: References

    - https://en.wikipedia.org/wiki/Comma-separated_values


.. _format_stockholm:

STOCKHOLM
---------

:Format: human readable
:Status: included
:Type: multiple sequence alignment

Stockholm format is a multiple sequence alignment format used by **Pfam** and **Rfam**
to store protein and RNA sequence alignments.

Here is a simple example::

    # STOCKHOLM 1.0
    #=GF ID    UPSK
    #=GF SE    Predicted; Infernal
    #=GF SS    Published; PMID 9223489
    #=GF RN    [1]
    #=GF RM    9223489
    #=GF RT    The role of the pseudoknot at the 3' end of turnip yellow mosaic
    #=GF RT    virus RNA in minus-strand synthesis by the viral RNA-dependent RNA
    #=GF RT    polymerase.
    #=GF RA    Deiman BA, Kortlever RM, Pleij CW;
    #=GF RL    J Virol 1997;71:5990-5996.

    AF035635.1/619-641             UGAGUUCUCGAUCUCUAAAAUCG
    M24804.1/82-104                UGAGUUCUCUAUCUCUAAAAUCG
    J04373.1/6212-6234             UAAGUUCUCGAUCUUUAAAAUCG
    M24803.1/1-23                  UAAGUUCUCGAUCUCUAAAAUCG
    #=GC SS_cons                   .AAA....<<<<aaa....>>>>
    //

A minimal well-formed Stockholm file should contain a header which states the
format and version identifier, currently '# STOCKHOLM 1.0', followed by the
sequences and corresponding unique sequence names::

    <seqname> <aligned sequence>
    <seqname> <aligned sequence>
    <seqname> <aligned sequence>

Mark-up lines may include any characters except whitespace. Use underscore ("_")
instead of space.::

    #=GF <feature> <Generic per-File annotation, free text>
    #=GC <feature> <Generic per-Column annotation, exactly 1 char per column>
    #=GS <seqname> <feature> <Generic per-Sequence annotation, free text>
    #=GR <seqname> <feature> <Generic per-Residue annotation, exactly 1 char per residue>

.. admonition:: Bioconvert conversions:

    :class:`~bioconvert.stockholm2clustal.STOCKHOLM2CLUSTAL`,
    :class:`~bioconvert.stockholm2phylip.STOCKHOLM2PHYLIP`


.. admonition:: References

    - https://en.wikipedia.org/wiki/Stockholm_format
    - http://scikit-bio.org/docs/0.5.0/generated/skbio.io.format.stockholm.html
    - pfam: https://en.wikipedia.org/wiki/Pfam



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

    - :class:`~bioconvert.bcf2vcf.BCF2VCF`
    - :class:`~bioconvert.bcf2wiggle.BCF2WIGGLE`
    - :class:`~bioconvert.vcf2bcf.VCF2BCF`
    - :class:`~bioconvert.vcf2bed.VCF2BED`
    - :class:`~bioconvert.vcf2wiggle.VCF2WIGGLE`
    - :class:`~bioconvert.vcf2plink.VCF2PLINK`
    - :class:`~bioconvert.vcf2bplink.VCF2BPLINK`


.. _format_wig:

WIG
---

See :ref:`format_wiggle`.

.. _format_wiggle:

WIGGLE (WIG)
------------

:Format: human readable
:Status: included
:Type: database-style


The wiggle (WIG) format is a format used for display of dense, continuous data such as GC percent.
Wiggle data elements must be equally sized.

Similar format such as the bedGraph format is also an older format used to display sparse data
or data that contains elements of varying size.

For speed and efficiency, wiggle data is usually stored in BIGWIG format.

Wiggle format is line-oriented. It is composed of declaration lines and data
lines. There are two options: **variableStep** and **fixedStep**.

The VariableStep format is used for data with irregular intervals between new data points,
and is the more commonly used wiggle format. The variableStep begins with a
declaration line and is followed by two columns containing chromosome positions
and data values::

    variableStep  chrom=chrN
    [span=windowSize]
      chromStartA  dataValueA
      chromStartB  dataValueB
      ... etc ...  ... etc ...

The declaration line starts with the word variableStep and is followed by a
specification for a chromosome. The optional span parameter (default: span=1)
allows data composed of contiguous runs of bases with the same data value to be
specified more succinctly. The span begins at each chromosome position specified
and indicates the number of bases that data value should cover. For example,
this variableStep specification::

    variableStep chrom=chr2
    300701 12.5
    300702 12.5
    300703 12.5
    300704 12.5
    300705 12.5

is equivalent to::

    variableStep chrom=chr2 span=5
    300701 12.5

The variableStep format becomes very inefficient when there are only a few data points per
1024 bases. If variableStep data points (i.e., chromStarts) are greater than
about 100 bases apart, it is advisable to use BedGraph format.

The **fixedStep** format is used for data with regular intervals between new data values and
is the more compact wiggle format. The fixedStep begins with a declaration line
and is followed by a single column of data values::

    fixedStep  chrom=chrN
    start=position  step=stepInterval
    [span=windowSize]
      dataValue1
      dataValue2
      ... etc ...

The declaration line starts with the word *fixedStep* and includes specifications
for chromosome, start coordinate, and step size. The span specification has the
same meaning as in variableStep format. For example, this fixedStep
specification::

    fixedStep chrom=chr3 start=400601 step=100
    11
    22
    33

displays the values 11, 22, and 33 as single-base regions on chromosome 3 at
positions 400601, 400701, and 400801, respectively. Adding span=5 to the
declaration line::

    fixedStep chrom=chr3 start=400601 step=100 span=5
    11
    22
    33

causes the values 11, 22, and 33 to be displayed as 5-base regions on chromosome
3 at positions 400601-400605, 400701-400705, and 400801-400805, respectively.

Note that for both variableStep and fixedStep formats, the same span must be
used throughout the dataset. If no span is specified, the default span of 1 is
used. As the name suggests, fixedStep wiggles require the same size step
throughout the dataset. If not specified, a step size of 1 is used.

Data values can be integer or real, postive or negative values.
Positions specified in the input data must be in numerical order.

.. warning::  BigWig files created from bedGraph format use "0-start, half-open"
    coordinates, but bigWigs that represent variableStep and fixedStep data are
    generated from wiggle files that use *1-start, fully-closed* coordinates. For
    example, for a chromosome of length N, the first position is 1 and the last
    position is N. For more information, see:


.. admonition:: Bioconvert conversions

    :class:`~bioconvert.wig2bed.WIG2BED`

.. admonition:: Reference

    - http://genome.ucsc.edu/goldenPath/help/wiggle.html

.. _format_xls:

XLS
---

:Format: human readable
:Type: database
:Status: included

Spreadsheet file format (Microsoft Excel file format).

Until 2007, Microsoft Excel used a proprietary binary file format
called Excel Binary File Format (.XLS). In Excel 2007, the Office Open XML was
introduced. We support the later formnat only.

With bioconvert you can convert an :ref:`format_xls` file into :ref:`format_csv` or :ref:`format_tsv` format. If several
sheets are to be found, you can select one or the other.


.. admonition:: Bioconvert conversions:

    :class:`~bioconvert.xls2csv.XLS2CSV`,
    :class:`~bioconvert.xlsx2csv.XLSZ2CSV`,

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

    :class:`~bioconvert.xls2csv.XLS2CSV`,
    :class:`~bioconvert.xlsx2csv.XLSX2CSV`

.. seealso::  :ref:`format_xls` format.

.. admonition::  References

    - https://en.wikipedia.org/wiki/Office_Open_XML



XMFA
----

:Format: human-readable
:Status: included
:Type: alignment

XMFA stands for eXtended Multi-FastA file format. The .alignment file contains
the complete genome alignment. This standard file format is also
used by other genome alignment systems that align sequences with rearrangements.

The XMFA file format supports the storage of several collinear sub-alignments,
each separated with an = sign, that constitute a single genome alignment. Each
sub-alignment consists of one FastA format sequence entry per genome where the
entry’s defline gives the strand (orientation) and location in the genome of the
sequence in the alignment.

Example (from darlinglab.org/mauve )::

    >seq_num:start1-end1 ± comments (sequence name, etc.)
    AC-TG-NAC--TG
    AC-TG-NACTGTG
    ...

    > seq_num:startN-endN ± comments (sequence name, etc.)
    AC-TG-NAC--TG
    AC-TG-NACTGTG
    ...
    = comments, and optional field-value pairs, i.e. score=12345

    > seq_num:start1-end1 ± comments (sequence name, etc.)
    AC-TG-NAC--TG
    AC-TG-NACTGTG
    ...

    > seq_num:startN-endN ± comments (sequence name, etc.)
    AC-TG-NAC--TG
    AC-TG-NACTGTG
    ...
    = comments, and optional field-value pairs, i.e. score=12345




.. admonition:: Bioconvert conversions

    :class:`~bioconvert.xmfa2phylip.XMFA2PHYLIP`


.. admonition:: References

    - http://darlinglab.org/mauve/user-guide/files.html


.. _format_yaml:

YAML
----

:Format: human-readable
:Status: included
:Type: database

YAML ("YAML Ain't Markup Language") is a human-readable data-serialization
language. It is commonly used for configuration files, but could be used in many
applications where data is being stored.


The full syntax cannot be described here. The full specification are available at the official site (https://yaml.org/refcard.html)

In brief:
- whitespace indentation is used to denote srtucture. Tab spaces are not allowed.
- **Comments** begin with the number sign #. Can start anywhere on a line.
- **List** are denoted by the - character with one member per line, or, enclosed in square brackets [ ] .
- **associated arrays** are represented with the colon space `:` in the form of *key:value*
- **strings** can be unquoted or quoted.

Example::

    # example of a yaml file
    - {name: Jean, age: 33}
    - name: Marie
      age : 32

    men:
        - Pierre
        - Jean
    women:
        - Marie


.. admonition:: Bioconvert conversions

    :class:`~bioconvert.json2yaml.JSON2YAML`,
    :class:`~bioconvert.yaml2json.JSON2YAML`.


.. admonition:: References

    - https://en.wikipedia.org/wiki/YAML
    - https://yaml.org/refcard.html

Others
------


ACE
~~~-

Human-readable file format used by the AceDB database, which is a genome database designed for
the handling of bioinformatics data. The data looks like::

    DNA : "HSFAU"
    ttccttccagctactgttccttccagc
    tactg

This format is obsolet and will not be included in **Bioconvert** for now.
BioPython seems to handle this format.

ASN1
~~~~

ASN.1 Abstract Syntax Notation One, is an International Standards
Organization (ISO) data representation format used to achieve interoperability.
It is formal notation used for describing data transmitted by
telecommunications protocols, regardless of language implementation and physical
representation of these data, whatever the application, whether complex or very
simple. NCBI uses ASN.1 for the storage and retrieval of data such as
nucleotide and protein sequences, structures, genomes, PubMed records, and more.

.. admonition:: Reference

    - https://www.ncbi.nlm.nih.gov/Structure/asn1.html
    - https://www.ncbi.nlm.nih.gov/IEB/ToolBox/SDKDOCS/ASNLIB.HTML#_AsnTool


.. _format_gcg:

GCG
~~~

:Format: human-readable
:Status: not included
:Type: sequence

GCG format contains exactly one sequence. It begins with
annotation lines and the start of the sequence is marked by a line ending with
two dot ("..") characters. This line also contains the sequence identifier, the
sequence length and a checksum. This format should only be used if the file was
created with the GCG package.

An example sequence in GCG format is::

    ID   AB000263 standard; RNA; PRI; 368 BP.
    XX
    AC   AB000263;
    XX
    DE   Homo sapiens mRNA for prepro cortistatin like peptide, complete cds.
    XX
    SQ   Sequence 368 BP;
    AB000263  Length: 368  Check: 4514  ..
           1  acaagatgcc attgtccccc ggcctcctgc tgctgctgct ctccggggcc acggccaccg
          61  ctgccctgcc cctggagggt ggccccaccg gccgagacag cgagcatatg caggaagcgg
         121  caggaataag gaaaagcagc ctcctgactt tcctcgcttg gtggtttgag tggacctccc
         181  aggccagtgc cgggcccctc ataggagagg aagctcggga ggtggccagg cggcaggaag
         241  gcgcaccccc ccagcaatcc gcgcgccggg acagaatgcc ctgcaggaac ttcttctgga
         301  agaccttctc ctcctgcaaa taaaacctca cccatgaatg ctcacgcaag tttaattaca
         361  gacctgaa

GVF
~~~

:Format: human-readable
:Status: not included
:Type: variant

The Genome Variation Format (GVF) is a very simple file format for describing
sequence_alteration features at nucleotide resolution relative to a reference
genome.

Example::

    ##gvf-version 1.10
    ##genome-build NCBI B36.3
    ##sequence-region chr16 1 88827254

    chr16 samtools SNV 49291141 49291141 . + . ID=ID_1;Variant_seq=A,G;Reference_seq=G;
    chr16 samtools SNV 49291360 49291360 . + . ID=ID_2;Variant_seq=G;Reference_seq=C;
    chr16 samtools SNV 49302125 49302125 . + . ID=ID_3;Variant_seq=T,C;Reference_seq=C;
    chr16 samtools SNV 49302365 49302365 . + . ID=ID_4;Variant_seq=G,C;Reference_seq=C;
    chr16 samtools SNV 49302700 49302700 . + . ID=ID_5;Variant_seq=T;Reference_seq=C;
    chr16 samtools SNV 49303084 49303084 . + . ID=ID_6;Variant_seq=G,T;Reference_seq=T;
    chr16 samtools SNV 49303156 49303156 . + . ID=ID_7;Variant_seq=T,C;Reference_seq=C;
    chr16 samtools SNV 49303427 49303427 . + . ID=ID_8;Variant_seq=T,C;Reference_seq=C;
    chr16 samtools SNV 49303596 49303596 . + . ID=ID_9;Variant_seq=T,C;Reference_seq=C;


.. admonition:: References

    - https://github.com/The-Sequence-Ontology/Specifications/blob/master/gvf.md


IG
--


The IntelliGenetics (IG) format is a sequence format. It can contain
several sequences, each consisting of a
number of comment lines that must begin with a semicolon (";"), a line with the
sequence name (it may not contain spaces!) and the sequence itself terminated
with the termination character '1' for linear or '2' for circular sequences.

An example sequence in IG format is::

    ; comment
    ; comment
    AB000263
    ACAAGATGCCATTGTCCCCCGGCCTCCTGCTGCTGCTGCTCTCCGGGGCCACGGCCACCGCTGCCCTGCC
    CCTGGAGGGTGGCCCCACCGGCCGAGACAGCGAGCATATGCAGGAAGCGGCAGGAATAAGGAAAAGCAGC
    CTCCTGACTTTCCTCGCTTGGTGGTTTGAGTGGACCTCCCAGGCCAGTGCCGGGCCCCTCATAGGAGAGG
    AAGCTCGGGAGGTGGCCAGGCGGCAGGAAGGCGCACCCCCCCAGCAATCCGCGCGCCGGGACAGAATGCC
    CTGCAGGAACTTCTTCTGGAAGACCTTCTCCTCCTGCAAATAAAACCTCACCCATGAATGCTCACGCAAG
    TTTAATTACAGACCTGAA1

PIR
~~~

:Format: human-readable
:Status: not included
:Type: variant

The PIR (Protein Informatics Resource) may contain contain several sequences.
A sequence in PIR format consists of One line starting with ">" character
followed by a 2-letter code describing the sequence type (P1, F1, DL, DC, RL,
RC, or XX), followed by a semicolon, followed by the sequence identification
code (the database ID-code). Then, one line containing a textual
description of the sequence and finally one or more lines containing the
sequence itself. The end of the sequence is marked by a "*"  character.

The PIR format is also often referred to as the NBRF format.

Example::

    >P1;CRAB_ANAPL
    Example protein sequence. Note the final * chraacter
    MDITIHNPLI RRPLFSWLAP SRIFDQIFGE HLQESELLPA SPSLSPFLMR
    SPIFRMPSWL ETGLSEMRLE KDKFSVNLDV KHFSPEELKV KVLGDMVEIH
    GKHEERQDEH GFIAREFNRK YRIPADVDPL TITSSLSLDG VLTVSAPRKQ
    SDVPERSIPI TREEKPAIAG AQRK*


.. admonition:: References

    - https://en.wikipedia.org/wiki/Protein_Information_Resource


- imgt    Unspecified (`*.txt`) This refers to the IMGT variant of the EMBL plain
  text file format.
- phd     PHD files are output from PHRED, used by PHRAP and CONSED for input.
- seqxml  Simple sequence XML file format.
- sff  Standard Flowgram Format (SFF) files produced by 454 sequencing.
  binary files produced by Roche 454 and IonTorrent/IonProton sequencing machines.
- swiss   Swiss-Prot aka UniProt format.
- uniprot-xml     UniProt XML format, successor to the plain text Swiss-Prot
  format.
- pdb2gmx: This program reads a .pdb (or .gro) file, reads some database files, adds
  hydrogens to the molecules and generates coordinates in GROMACS (GROMOS), or
  optionally .pdb, format and a topology in GROMACS format. See
  http://manual.gromacs.org/archive/4.6.7/online/pdb2gmx.html for details.
  this tool is already quite complete and will not be provided for now.
- rfam: https://en.wikipedia.org/wiki/Rfam

