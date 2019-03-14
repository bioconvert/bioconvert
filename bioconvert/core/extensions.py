###########################################################################
# Bioconvert is a project to facilitate the interconversion               #
# of life science data from one format to another.                        #
#                                                                         #
# Authors: see CONTRIBUTORS.rst                                           #
# Copyright © 2018  Institut Pasteur, Paris and CNRS.                     #
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

from easydev import AttrDict


extensions = {
    'abi': ["abi", "ab1"],                      # sequence
    'bam': ["bam"],                             # alignment
    'bcf': ["bcf"],                             # variant
    'bed': ["bed"],                             # misc
    'bedgraph': ["bedgraph"],
    'bigwig': ["bigwig"],
    'bz2': ['bz2'],
    'cdao': ["cdao"],                           # phylo
    'cram': ["cram"],                           # alignment
    'clustal': ["clustal", "aln", "clw"],       # phylo
    'csv': ["csv"],
    'dsrc': ['dsrc'],
    'embl': ['embl'],                           # annotation/sequence
    'fasta': ["fasta", "fa", "fst"],            # sequence
    'fastq': ["fastq", "fq"],                   # sequence
    'genbank': ['genbank', 'gbk', "gb"],        # annotation/sequence
    'gfa': ['gfa'],                             # assembly
    'gff': ['gff3', 'gff'],                     # gff3 so far
    'gz': ['gz'],
    'json': ['json'],                           # misc
    'newick': ["newick", "nw", "nhx", "nwk"],   # phylo
    'nexus': ["nexus", "nx", "nex", "nxs"],     # phylo
    'paf': ['paf'],                             # assembly
    'phylip': ['phy', 'ph', 'phylip'],          # phylo
    'phyloxml': ['phyloxml', 'xml'],            # phylo
    'sam': ["sam"],                             # alignement
    'sra': ["sra"],                             # sra format
    'stockholm': ['sto', 'sth', 'stockholm'],   # alignment
    'vcf': ['vcf'],                             # variant
    'twobit': ['2bit'],                         # sequence
    'tsv': ["tsv"],
    'yaml': ['yaml', 'YAML'],                   # misc
    'maf': ["maf"]                              # !! this is MIRA format, not mutation alignment format
}

extensions = AttrDict(**extensions)


# nexml   *.xml
# phyloxml    *.xml

"""
ace     *.ace   1.47    No  1.52    Reads the contig sequences from an ACE assembly file. Uses Bio.Sequencing.Ace internally   
clustal     *.aln   1.43    1.43    No  The alignment format of Clustal X and
Clustal W.  CLUSTAL format is recognised by the word CLUSTAL at the beginning of
the file.

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

ig  Unspecified (*.txt)     1.47    No  1.52    This refers to the
IntelliGenetics file format, apparently the same as the MASE alignment format.

imgt    Unspecified (*.txt)     1.56    1.56    1.56    This refers to the IMGT
variant of the EMBL plain text file format.


phd     *.phd   1.46    1.52    1.52    PHD files are output from PHRED,
used by PHRAP and CONSED for input. Uses Bio.Sequencing.Phd internally.

pir     *.pir   1.48    No  1.52    A "FASTA like" format introduced by the
National Biomedical Research Foundation (NBRF) for the Protein Information
Resource (PIR) database, now part of UniProt.

sff     *.sff   1.54    1.54    1.54    Standard Flowgram Format (SFF) binary
files produced by Roche 454 and IonTorrent/IonProton sequencing machines.

swiss   *.sw    1.43    No  1.52    Swiss-Prot aka UniProt format. Uses
Bio.SwissProt internally. See also the UniProt XML format.

qual    *.qual  1.50    1.50    1.52    Qual files are a bit like FASTA files
but instead of the sequence, record space separated integer sequencing values
as PHRED quality scores. A matched pair of FASTA and QUAL files are often used
as an alternative to a single FASTQ file.

uniprot-xml     *.xml   1.56    No  1.56    UniProt XML format, successor to
the plain text Swiss-Prot format.

"""
