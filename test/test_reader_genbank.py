from bioconvert import bioconvert_data
from bioconvert.io.genbank import Genbank

import pytest


def test_load():
    infile = bioconvert_data("test_genbank.gbk")
    reader = Genbank(infile)

    for gbk in reader.read():
        assert len(gbk) == 9

        # LOCUS tests
        loc = gbk["LOCUS"]
        assert len(loc) == 5

        assert loc["id"] == "AF068625"
        assert loc["length"] == "200"
        assert loc["mol_type"] == "mRNA linear"
        assert loc["genbank_div"] == "ROD"
        assert loc["date"] == "06-DEC-1999"

        # DEFINITION
        assert gbk["DEFINITION"] == "Mus musculus DNA cytosine-5 methyltransferase 3A (Dnmt3a) mRNA, complete cds."

        # ACCESSION
        # 'ACCESSION': {'id': 'AF068625', 'other': 'REGION: 1..200'}
        acs = gbk["ACCESSION"]
        assert len(acs) == 2
        assert acs["id"] == "AF068625"
        assert acs["other"] == "REGION: 1..200"

        # VERSION
        # VERSION     AF068625.2  GI:6449467
        # 'id': 'AF068625.2', 'GI': '6449467'
        ver = gbk["VERSION"]
        assert len(ver) == 2
        assert ver["id"] == "AF068625.2"
        assert ver["GI"] == "6449467"

        # SOURCE
        # 'SOURCE': {'short': 'Mus musculus (house mouse)',
        src = gbk["SOURCE"]
        assert len(src) == 2
        assert src["short"] == "Mus musculus (house mouse)"

        # SOURCE ORGANISM
        orga = src["ORGANISM"]
        assert len(orga) == 2
        assert orga["name"] == "Mus musculus" 
        assert orga["lineage"] == "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Glires; Rodentia; Sciurognathi; Muroidea; Muridae; Murinae; Mus."

        # 'REFERENCE': [, , ],
        ref = gbk["REFERENCE"]
        assert len(ref) == 3
        assert len(ref[0]) == 4
        assert len(ref[1]) == 3
        assert len(ref[2]) == 4

        # {'AUTHORS': 'Okano,M., Xie,S. and Li,E.', 'TITLE': 'Cloning and characterization of a family of novel mammalian DNA', 'JOURNAL': 'Nat. Genet. 19 (3), 219-220 (1998)', 'PUBMED': '9662389'}
        assert ref[0]["AUTHORS"] == "Okano,M., Xie,S. and Li,E."
        assert ref[0]["TITLE"] == "Cloning and characterization of a family of novel mammalian DNA (cytosine-5) methyltransferases"
        assert ref[0]["JOURNAL"] == "Nat. Genet. 19 (3), 219-220 (1998)"
        assert ref[0]["PUBMED"] == "9662389"

        # {'AUTHORS': 'Xie,S., Okano,M. and Li,E.', 'TITLE': 'Direct Submission', 'JOURNAL': 'Submitted (28-MAY-1998) CVRC, Mass. Gen. Hospital, 149 13th Street, Charlestown, MA 02129, USA'}
        assert ref[1]["AUTHORS"] == "Xie,S., Okano,M. and Li,E."
        assert ref[1]["TITLE"] == "Direct Submission"
        assert ref[1]["JOURNAL"] == "Submitted (28-MAY-1998) CVRC, Mass. Gen. Hospital, 149 13th Street, Charlestown, MA 02129, USA"

        # {'AUTHORS': 'Okano,M., Chijiwa,T., Sasaki,H. and Li,E.', 'TITLE': 'Direct Submission', 'JOURNAL': 'Submitted (04-NOV-1999) CVRC, Mass. Gen. Hospital, 149 13th Street,', 'REMARK': 'Sequence update by submitter'}
        assert ref[2]["AUTHORS"] == "Okano,M., Chijiwa,T., Sasaki,H. and Li,E."
        assert ref[2]["TITLE"] == "Direct Submission"
        assert ref[2]["JOURNAL"] == "Submitted (04-NOV-1999) CVRC, Mass. Gen. Hospital, 149 13th Street, Charlestown, MA 02129, USA"
        assert ref[2]["REMARK"] == "Sequence update by submitter"

        # 'COMMENT': 'On Nov 18, 1999 this sequence version replaced gi:3327977.',
        assert gbk["COMMENT"] == "On Nov 18, 1999 this sequence version replaced gi:3327977."

        # 'FEATURES': [, ], 
        fea = gbk["FEATURES"]
        assert len(fea) == 2

        # {'type': 'source', 'position': '1..200', 'organism': 'Mus musculus', 'mol_type': 'mRNA', 'db_xref': 'taxon:10090', 'chromosome': '12', 'map': '4.0 cM'}
        assert fea[0]["type"] == "source"
        assert fea[0]["position"] == "1..200"
        assert fea[0]["organism"] == "Mus musculus"
        assert fea[0]["mol_type"] == "mRNA"
        assert fea[0]["db_xref"] == "taxon:10090"
        assert fea[0]["chromosome"] == "12"
        assert fea[0]["map"] == "4.0 cM"

        # {'type': 'gene', 'position': '1..>200', 'gene': 'Dnmt3a'}
        assert fea[1]["type"] == "gene"
        assert fea[1]["position"] == "1..>200"
        assert fea[1]["gene"] == "Dnmt3a"

        # 'ORIGIN': 'gaattccggcctgctgccgggccgcccgacccgccgggccacacggcagagccgcctgaagcccagcgctgaggctgcacttttccgagggcttgacatcagggtctatgtttaagtcttagctcttgcttacaaagaccacggcaattccttctctgaagccctcgcagccccacagcgccctcgcagccccagcctgc'}
        assert len(gbk["ORIGIN"]) == int(gbk["LOCUS"]["length"]) == 200
        assert gbk["ORIGIN"] == "gaattccggcctgctgccgggccgcccgacccgccgggccacacggcagagccgcctgaagcccagcgctgaggctgcacttttccgagggcttgacatcagggtctatgtttaagtcttagctcttgcttacaaagaccacggcaattccttctctgaagccctcgcagccccacagcgccctcgcagccccagcctgc"

