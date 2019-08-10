from bioconvert import bioconvert_data
from bioconvert.io.gff3 import Gff3

import pytest


expected_values = [{'seqid': 'ctg123', 'type': 'gene', 'start': 1000, 'stop': 9000, 'strand': '+', 'attributes': {'ID': 'gene00001', 'Name': 'EDEN'}}
, {'seqid': 'ctg123', 'type': 'TF_binding_site', 'start': 1000, 'stop': 1012, 'strand': '+', 'attributes': {'ID': 'tfbs00001', 'Parent': 'gene00001'}}
, {'seqid': 'ctg123', 'type': 'mRNA', 'start': 1050, 'stop': 9000, 'strand': '+', 'attributes': {'ID': 'mRNA00001', 'Parent': 'gene00001', 'Name': 'EDEN.1'}}
, {'seqid': 'ctg123', 'type': 'mRNA', 'start': 1050, 'stop': 9000, 'strand': '+', 'attributes': {'ID': 'mRNA00002', 'Parent': 'gene00001', 'Name': 'EDEN.2'}}
, {'seqid': 'ctg123', 'type': 'mRNA', 'start': 1300, 'stop': 9000, 'strand': '+', 'attributes': {'ID': 'mRNA00003', 'Parent': 'gene00001', 'Name': 'EDEN.3'}}
, {'seqid': 'ctg123', 'type': 'exon', 'start': 1300, 'stop': 1500, 'strand': '+', 'attributes': {'ID': 'exon00001', 'Parent': 'mRNA00003'}}
, {'seqid': 'ctg123', 'type': 'exon', 'start': 1050, 'stop': 1500, 'strand': '+', 'attributes': {'ID': 'exon00002', 'Parent': 'mRNA00001,mRNA00002'}}
, {'seqid': 'ctg123', 'type': 'exon', 'start': 3000, 'stop': 3902, 'strand': '+', 'attributes': {'ID': 'exon00003', 'Parent': 'mRNA00001,mRNA00003'}}
, {'seqid': 'ctg123', 'type': 'exon', 'start': 5000, 'stop': 5500, 'strand': '+', 'attributes': {'ID': 'exon00004', 'Parent': 'mRNA00001,mRNA00002,mRNA00003'}}
, {'seqid': 'ctg123', 'type': 'exon', 'start': 7000, 'stop': 9000, 'strand': '+', 'attributes': {'ID': 'exon00005', 'Parent': 'mRNA00001,mRNA00002,mRNA00003'}}
, {'seqid': 'ctg123', 'type': 'CDS', 'start': 1201, 'stop': 1500, 'strand': '+', 'phase': 0, 'attributes': {'ID': 'cds00001', 'Parent': 'mRNA00001', 'Name': 'edenprotein.1'}}
, {'seqid': 'ctg123', 'type': 'CDS', 'start': 3000, 'stop': 3902, 'strand': '+', 'phase': 0, 'attributes': {'ID': 'cds00001', 'Parent': 'mRNA00001', 'Name': 'edenprotein.1'}}
, {'seqid': 'ctg123', 'type': 'CDS', 'start': 5000, 'stop': 5500, 'strand': '+', 'phase': 0, 'attributes': {'ID': 'cds00001', 'Parent': 'mRNA00001', 'Name': 'edenprotein.1'}}
, {'seqid': 'ctg123', 'type': 'CDS', 'start': 7000, 'stop': 7600, 'strand': '+', 'phase': 0, 'attributes': {'ID': 'cds00001', 'Parent': 'mRNA00001', 'Name': 'edenprotein.1'}}
, {'seqid': 'ctg123', 'type': 'CDS', 'start': 1201, 'stop': 1500, 'strand': '+', 'phase': 0, 'attributes': {'ID': 'cds00002', 'Parent': 'mRNA00002', 'Name': 'edenprotein.2'}}
, {'seqid': 'ctg123', 'type': 'CDS', 'start': 5000, 'stop': 5500, 'strand': '+', 'phase': 0, 'attributes': {'ID': 'cds00002', 'Parent': 'mRNA00002', 'Name': 'edenprotein.2'}}
, {'seqid': 'ctg123', 'type': 'CDS', 'start': 7000, 'stop': 7600, 'strand': '+', 'phase': 0, 'attributes': {'ID': 'cds00002', 'Parent': 'mRNA00002', 'Name': 'edenprotein.2'}}
, {'seqid': 'ctg123', 'type': 'CDS', 'start': 3301, 'stop': 3902, 'strand': '+', 'phase': 0, 'attributes': {'ID': 'cds00003', 'Parent': 'mRNA00003', 'Name': 'edenprotein.3'}}
, {'seqid': 'ctg123', 'type': 'CDS', 'start': 5000, 'stop': 5500, 'strand': '+', 'phase': 1, 'attributes': {'ID': 'cds00003', 'Parent': 'mRNA00003', 'Name': 'edenprotein.3'}}
, {'seqid': 'ctg123', 'type': 'CDS', 'start': 7000, 'stop': 7600, 'strand': '+', 'phase': 1, 'attributes': {'ID': 'cds00003', 'Parent': 'mRNA00003', 'Name': 'edenprotein.3'}}
, {'seqid': 'ctg123', 'type': 'CDS', 'start': 3391, 'stop': 3902, 'strand': '+', 'phase': 0, 'attributes': {'ID': 'cds00004', 'Parent': 'mRNA00003', 'Name': 'edenprotein.4'}}
, {'seqid': 'ctg123', 'type': 'CDS', 'start': 5000, 'stop': 5500, 'strand': '+', 'phase': 1, 'attributes': {'ID': 'cds00004', 'Parent': 'mRNA00003', 'Name': 'edenprotein.4'}}
, {'seqid': 'ctg123', 'type': 'CDS', 'start': 7000, 'stop': 7600, 'strand': '+', 'phase': 1, 'attributes': {'ID': 'cds00004', 'Parent': 'mRNA00003', 'Name': 'edenprotein.4'}}]


def test_load():
    infile = bioconvert_data("GFF3/gff3_example.gff")
    reader_gff3 = Gff3(infile)

    for expected_value, file_value in zip(expected_values, reader_gff3.read()):
        assert expected_value == file_value

