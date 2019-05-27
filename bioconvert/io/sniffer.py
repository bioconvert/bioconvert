from bioconvert.core.extensions import extensions



env =    PYTHONHASHSEED=0

class Sniffer(object):
    """

    :: 

        >>> from bioconvert import Sniffer
        >>> s =  Sniffer()
        >>> s.sniff("test.clustal")
        "clustal"
        

    """
    formats = sorted(extensions.keys())

    def __init__(self):
        pass

    def sniff(self, filename):
        """Return first frmt found to be compatible with the input file"""
        for frmt in self.formats:
            func = getattr(self, "is_{}".format(frmt))
            try:
                ret = func(filename)
                if ret:
                    return frmt
            except NotImplementedError:
                pass
            #except Exception as err:
            #    raise(err)
        return None

    def _is_blank_line(self, line):
        line = line.strip()
        if len(line) == 0:
            return True
        else:
            return False

    def is_abi(self, filename):
        raise NotImplementedError

    def is_bam(self, filename):
        raise NotImplementedError

    def is_bcf(self, filename):
        raise NotImplementedError

    def is_bed(self, filename):
        raise NotImplementedError

    def is_bedgraph(self, filename):
        raise NotImplementedError

    def is_bigwig(self, filename):
        raise NotImplementedError

    def is_bigbed(self, filename):
        raise NotImplementedError

    def is_bplink(self, filename):
        raise NotImplementedError

    def is_bz2(self, filename):
        raise NotImplementedError

    def is_cdao(self, filename):
        raise NotImplementedError

    def is_csv(self, filename):
        raise NotImplementedError

    def is_cram(self, filename):
        raise NotImplementedError

    def is_clustal(self, filename):
        with open(filename, "r") as fin:
            try:
                line = fin.readline().strip()
                if self._is_blank_line(line):
                    pass
                elif line.startswith("CLUSTAL"):
                    return True
            except:
                return False

    def is_dsrc(self, filename):
        raise NotImplementedError

    def is_embl(self, filename):
        raise NotImplementedError

    def is_fasta(self, filename):
        raise NotImplementedError

    def is_fastq(self, filename):
        raise NotImplementedError

    def is_genbank(self, filename):

        with open(filename, "r") as fin:
            try:
                line = fh.readline().strip()
                if self._is_blank_line(line):
                    pass
                else:
                    data = line.split()
                    data[0] in ['LOCUS']
                return True
            except:
                return False

    def is_gfa(self, filename):
        raise NotImplementedError

    def is_gff2(self, filename):
        raise NotImplementedError

    def is_gff3(self, filename):
        raise NotImplementedError

    def is_gz(self, filename):
        raise NotImplementedError

    def is_json(self, filename):
        raise NotImplementedError

    def is_maf(self, filename):
        raise NotImplementedError

    def is_newick(self, filename):
        raise NotImplementedError

    def is_nexus(self, filename):
        raise NotImplementedError

    def is_ods(self, filename):
        raise NotImplementedError

    def is_paf(self, filename):
        raise NotImplementedError

    def is_phylip(self, filename):
        raise NotImplementedError

    def is_phyloxml(self, filename):
        raise NotImplementedError

    def is_plink(self, filename):
        raise NotImplementedError

    def is_qual(self, filename):
        raise NotImplementedError

    def is_sam(self, filename):
        raise NotImplementedError

    def is_scf(self, filename):
        raise NotImplementedError

    def is_sra(self, filename):
        raise NotImplementedError

    def is_stockholm(self, filename):
        raise NotImplementedError

    def is_twobit(self, filename):
        raise NotImplementedError

    def is_tsv(self, filename):
        raise NotImplementedError

    def is_vcf(self, filename):
        raise NotImplementedError

    def is_wiggle(self, filename):
        raise NotImplementedError

    def is_wig(self, filename):
        raise NotImplementedError

    def is_xls(self, filename):
        raise NotImplementedError

    def is_xlsx(self, filename):
        raise NotImplementedError

    def is_xmfa(self, filename):
        raise NotImplementedError

    def is_yaml(self, filename):
        raise NotImplementedError


