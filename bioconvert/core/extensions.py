###########################################################################
# Bioconvert is a project to facilitate the interconversion               #
# of life science data from one format to another.                        #
#                                                                         #
# Copyright © 2018-2022  Institut Pasteur, Paris and CNRS.                #
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
#                                                                         #
# Repository: https://github.com/bioconvert/bioconvert                    #
# Documentation: http://bioconvert.readthedocs.io                         #
###########################################################################
"""List of formats and associated extensions"""


class AttrDict(dict):
    """Copy from easydev package."""

    def __init__(self, **kwargs):
        dict.__init__(self, kwargs)
        self.__dict__ = self
        self.update(kwargs)

    def update(self, content):
        """See class/constructor documentation for details

        :param dict content: a valid dictionary
        """
        # accepts dict and attrdict classes
        try:
            from collections import OrderedDict
        except:
            OrderedDict = AttrDict

        if content.__class__ not in [dict, OrderedDict, AttrDict]:
            raise TypeError

        for k, v in content.items():
            if v.__class__ not in [dict, AttrDict, OrderedDict]:
                # fixme copy ?
                self[k] = v
            else:
                self[k] = AttrDict(**v)


# Formats can be of type
# - sequence
# - alignment
# - binary
# - compression
# - database
# - variant

#: List of formats and their extensions included in Bioconvert
extensions = {
    "abi": ["abi", "ab1"],  # sequence
    "agp": ["agp"],  # assembly
    "bam": ["bam"],  # alignment
    "bcf": ["bcf"],  # variant
    "bed": ["bed"],  # database
    "bedgraph": ["bedgraph", "bg"],  # database
    "bigwig": ["bigwig", "bw"],  # database
    "bigbed": ["bb", "bigbed"],
    "bz2": ["bz2"],  # compression
    "bplink": ["bplink"],
    "cdao": ["cdao"],  # phylo
    "cram": ["cram"],  # alignment
    "clustal": ["clustal", "aln", "clw"],  # phylo
    "cov": ["cov"],  # coverage (chrom name,  pos, depth)
    "csv": ["csv"],  # database
    "dsrc": ["dsrc"],  # compression
    "embl": ["embl"],  # annotation/sequence
    "ena": ["ena"],
    "faa": ["faa", "mpfa", "aa"],  # fasta multiple amino acid
    "fast5": ["fast5"],
    "fasta": ["fasta", "fa", "fst"],  # sequence
    "fastq": ["fastq", "fq"],  # sequence
    "genbank": ["genbank", "gbk", "gb"],  # annotation/sequence
    "gfa": ["gfa"],  # assembly
    "gff2": ["gff"],
    "gff3": ["gff3"],  # annotation
    "gtf": ["gtf"],  # annotation
    "gz": ["gz"],
    "json": ["json"],  # database
    "maf": ["maf"],  # !! this is MIRA format, not mutation alignment format
    "newick": ["newick", "nw", "nhx", "nwk"],  # phylo
    "nexus": ["nexus", "nx", "nex", "nxs"],  # phylo
    "ods": ["ods"],  # database
    "paf": ["paf"],  # assembly
    "pdb": ["pdb"],
    "phylip": ["phy", "ph", "phylip"],  # phylo
    "phyloxml": ["phyloxml", "xml"],  # phylo
    "plink": ["plink"],
    "pod5": ["pod5"],
    "qual": ["qual"],  # seauence
    "sam": ["sam"],  # alignement
    "scf": ["scf"],  # alignement
    "sra": ["sra"],  # sra format
    "stockholm": ["sto", "sth", "stk", "stockholm"],  # alignment
    "twobit": ["2bit"],  # sequence
    "tsv": ["tsv"],  # database
    "vcf": ["vcf"],  # variant
    "wiggle": ["wig", "wiggle"],
    "wig": ["wig"],
    "xls": ["xls"],  # database
    "xlsx": ["xlsx"],  # database
    "xmfa": ["xmfa"],
    "yaml": ["yaml", "YAML"],  # database
}

extensions = AttrDict(**extensions)
