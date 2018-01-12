


from easydev import AttrDict


extensions = {
    'abi': ["abi", "ab1"],                      # sequence
    'cdao': ["cdao"],                           # phylo
    'clustal':["clustal", "aln", "clw"],        # phylo
    'fastq': ["fasta", "fq"],                   # sequence
    'newick': ["newick", "nw", "nhx", "nwk"],   # phylo
    'nexus': ["nexus", "nx", "nex", "nxs"],     # phylo
    'phylip': ['phy', 'ph', 'phylip'],          # phylo
    'stockholm': ['sto', 'sth', 'stockholm'],   # alignment

    }

extensions = AttrDict(**extensions)


# nexml   *.xml   
# phyloxml    *.xml  

"""
ace     *.ace   1.47    No  1.52    Reads the contig sequences from an ACE assembly file. Uses Bio.Sequencing.Ace internally   
clustal     *.aln   1.43    1.43    No  The alignment format of Clustal X and
Clustal W.  CLUSTAL format is recognised by the word CLUSTAL at the beginning of
the file.

embl    Unspecified (*.txt)     1.43    1.54    1.52    The EMBL flat file
format. Uses Bio.GenBank internally.    
fasta   *.fasta, *.fas, *.fa, *.seq, *.fsa, *.fna, *.ffn, *.faa, *.frn  1.43
1.43    1.52    This refers to the input FASTA file format introduced for Bill
Pearson's FASTA tool, where each record starts with a ">" line. Resulting
sequences have a generic alphabet by default.   There is no standard file
extension for a text file containing FASTA formatted sequences. Although their
is a plethora of ad-hoc file extensions. See this article for details and
explanation.

fastq-solexa    *.fq, *.fastq   1.50    1.50    1.52    FASTQ files are a bit
like FASTA files but also include sequencing qualities. In
Biopython,"fastq-solexa" refers to the original Solexa/Illumina style FASTQ
files which encode Solexa qualities using an ASCII offset of 64. See also what
we call the "fastq-illumina" format.    There is no standard file extension for
a FASTQ file, but .fq and .fastq, are commonly used. There are different FASTQ
formats for different sequencing technologies.

fastq-illumina  *.fq, *.fastq   1.51    1.51    1.52    FASTQ files are a bit
like FASTA files but also include sequencing qualities. In
Biopython,"fastq-illumina" refers to early Solexa/Illumina style FASTQ files
(from pipeline version 1.3 to 1.7) which encode PHRED qualities using an ASCII
offset of 64. For good quality reads, PHRED and Solexa scores are approximately
equal, so the "fastq-solexa" and "fastq-illumina" variants are almost
equivalent.     There is no standard file extension for a FASTQ file, but .fq
and .fastq, are commonly used. There are different FASTQ formats for different
sequencing technologies.

genbank or gb   *.gbk, *.gb, *.gpff     1.43    1.48 / 1.51     1.52    The
GenBank or GenPept flat file format. Uses Bio.GenBank internally for parsing.
Biopython 1.48 to 1.50 wrote basic GenBank files with only minimal annotation,
while 1.51 onwards will also write the features table (see Bug 2294).   See this
article.

ig  Unspecified (*.txt)     1.47    No  1.52    This refers to the
IntelliGenetics file format, apparently the same as the MASE alignment format.  

imgt    Unspecified (*.txt)     1.56    1.56    1.56    This refers to the IMGT
variant of the EMBL plain text file format.     


phd     *.phd   1.46    1.52    1.52    PHD files are output from PHRED, used by
PHRAP and CONSED for input. Uses Bio.Sequencing.Phd internally.     

pir     *.pir   1.48    No  1.52    A "FASTA like" format introduced by the
National Biomedical Research Foundation (NBRF) for the Protein Information
Resource (PIR) database, now part of UniProt.   

sff     *.sff   1.54    1.54    1.54    Standard Flowgram Format (SFF) binary
files produced by Roche 454 and IonTorrent/IonProton sequencing machines.   

swiss   *.sw    1.43    No  1.52    Swiss-Prot aka UniProt format. Uses
Bio.SwissProt internally. See also the UniProt XML format.  

qual    *.qual  1.50    1.50    1.52    Qual files are a bit like FASTA files
but instead of the sequence, record space separated integer sequencing values as
PHRED quality scores. A matched pair of FASTA and QUAL files are often used as
an alternative to a single FASTQ file.  

uniprot-xml     *.xml   1.56    No  1.56    UniProt XML format, successor to the
plain text Swiss-Prot format.

"""
