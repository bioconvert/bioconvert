
.. _glossary:

Glossary
========


Note that formats mentionned here below have dedicated description in the
:ref:`formats` section. 

.. glossary::
    :sorted:

    ABI

        File format produced by ABI sequencing machines. Contains the trace data
        which includes probabilities of the four nucleotides. See
        the :ref:`format_abi` format page for details.

    ASQG

       The ASQG format describes an assembly graph. Each line is a tab-delimited
       record. The first field in each record describes the record type. See the
       :ref:`format_asqg` page for details.

    BAI

        The index file related to file generated in the BAM format. (This is a
        non-standard file type.) See the :ref:`format_bai` page for details.

    BAM

        Binary version of the Sequence Alignment Map (SAM) format. See the
        :ref:`format_bam` format page for details. 

    BCF

        Binary version of the Variant Call Format (VCF).
        See :ref:`format_bcf` page for details.

    BCL

        BCL is the raw format used by Illumina sequencers. See the :ref:`format_bcl` format 
        page for details.

    BED

        BEDGRAPH/BED format is line-oriented and allows display of continuous-valued
        data. Similar to WIG format.
        See the :ref:`format_bed` format page for details.


    BED3

        Variants of the BED format with 4 columns storing the track name,
        start and end positions and values.
        See the :ref:`format_bed4` format page for details.

    BED4

        Variants of the BED format with 4 columns storing the track name,
        start and end positions and values.
        See the :ref:`format_bed4` format page for details.

    BEDGRAPH

        BEDGRAPH/BED format is line-oriented and allows display of continuous-valued
        data. Similar to WIG format.
        See the :ref:`format_bed` format page for details.

    BIGBED

        An indexed binary version of a BED file
        See :ref:`format_bigbed` page for details.

    BPLINK

        Binary version of the PlINK forat used for analyzing genotypic data 
        for Genome-wide Association Studies (GWAS). 
        See :ref:`format_plink_binary` page for details.

    BZ2

        **bzip2** is a file compression program that uses the Burrows–Wheeler algorithm. 
        Extension is usually .bz2 See :ref:`format_bz2` page for details.

    BIGWIG

        Indexed binary version of the Wiggle format.
        See :ref:`format_bigwig` page for details.

    CLUSTAL

        The alignment format of Clustal X and Clustal W. See
        :ref:`format_clustal` page for details.

    COV

        A bioconvert format to store coverage in the form of a 3 column 
        tab-tabulated file. See :ref:`format_cov` page for details.

    CRAM

        A more compact version of BAM files used to store Sequence Alignment 
        Map (SAM) format. See :ref:`format_cram` page for details.

    CSV

        A comma-separated values format is a delimited text file that uses a
        comma to separate values. See :ref:`format_csv` format page for
        details.

    DSRC

        A compression tool dedicated to FastQ files
        See :ref:`format_dsrc` page for details.

    EMBL

        EMBL Flat File Format.
        See :ref:`format_embl` page for details.

    FAA

        FASTA-formatted sequence files containing amino acid sequences
        See :ref:`format_faa` page for details.
 
    FASTA

        FASTA-formatted sequence files contain either nucleic acid sequence
        (such as DNA) or protein sequence information. FASTA files can also store multiple
        sequences in a single file. See :ref:`format_fasta` page for details.

    FASTQ

        FASTQ-formatted sequence files are used to represent high-throughput
        sequencing data, where each read is described by a name, its sequence,
        and its qualities. See :ref:`format_fastq` page for details.

    GFA

        Graphical Fragment Assembly format. https://github.com/GFA-spec/GFA-spec

    GFF2

        General Feature Format, used for describing genes and other features
        associated with DNA, RNA and Protein sequences.
        See :ref:`format_gff` page for details.

    GFF3

        General Feature Format, used for describing genes and other features
        associated with DNA, RNA and Protein sequences.
        http://genome.ucsc.edu/FAQ/FAQformat#format3
        See :ref:`format_gff` page for details.

    GENBANK

        GenBank Flat File Format.
        See :ref:`format_genbank` page for details.

    GZ

        **gzip** is a file compression program based on the DEFLATE algorithm. 
        See :ref:`format_gz` page for details.

    JSON

        A human-readable data serialization language commonly used in
        configuration files. See :ref:`format_json` page for details.

    MAF

        A human-readable multiple alignment format. 
        See :ref:`format_maf` page for details.

    NEXUS

        Plain text minimal format used to store multiple alignment and 
        phylogenetic trees. See :ref:`format_nexus` page for details.

    NEWICK

        Plain text minimal format used to store phylogenetic tree.
        See :ref:`format_newick` page for details.

    PAF

        PAF is a text format describing the approximate mapping positions
        between two set of sequences.

    PHYLIP

        Plain text format to store a multiple sequence alignment.
        See :ref:`format_phylip` page for details.

    PHYLOXML

        XML format to store a multiple sequence alignment.
        See :ref:`format_phyloxml` page for details.

    PLINK

        Format used for analyzing genotypic data for Genome-wide Association
        Studies (GWAS). See :ref:`format_plink_flat` page for details.

    QUAL

        Sequence of qualities associated with a sequence of nucleotides.
        Associated with FastA file, the original FastQ file can be built back.
        See :ref:`format_qual` page for details.

    SAM

        Sequence Alignment Map is a generic nucleotide alignment format that
        describes the alignment of query sequences or sequencing reads to a reference
        sequence or assembly. See :ref:`format_sam` page for details.

    SCF

        Standard Chromatogram Format, a binary
        chromatogram format described in Staden package documentation SCF file format.

    SRA

        The Sequence Read Archive (SRA) is a website that stores
        sequencing data at https://www.ncbi.nlm.nih.gov/sra
        It is not a format per se. See :ref:`format_sra` page for details.

    STOCKHOLM

        Stockholm format is a multiple sequence alignment format used to store 
        multiple sequence alignment. See :ref:`format_stockholm` page for details.

    TSV

        A tab-separated values format is a delimited text file that uses a
        tab character to separate values. See :ref:`format_tsv` format page for
        details.

    TWOBIT

        **2bit** file stores multiple DNA sequences (up to 4 Gb total) in a
        compact randomly-accessible format. The file contains masking information 
        as well as the DNA itself. See :ref:`format_twobit` format page for
        details.

    VCF

        Variant Call Format (VCF) is a flexible and extendable format for 
        storing variation in sequences such as single nucleotide variants,
        insertions/deletions, copy number variants and structural variants. 
        See :ref:`format_vcf` page for details.

    WIG

        Synonym for the wiggle (WIG) format. See :ref:`format_wig`.

    WIGGLE

        The wiggle (WIG) format stores dense, continuous data such as GC percent, 
        probability scores, and transcriptome data. See :ref:`format_wig` page
        for details.

    XLS

        Spreadsheet file format (Microsoft Excel file format). 
        See :ref:`format_xls` page for details.

    XLSX

        Spreadsheet file format defined in the Office Open XML specification.
        See :ref:`format_xlsx` page for details.


    XMFA

        TODO

    YAML

        A human-readable data serialization language commonly used in
        configuration files. See https://en.wikipedia.org/wiki/YAML
        See :ref:`format_yaml` page for details.


