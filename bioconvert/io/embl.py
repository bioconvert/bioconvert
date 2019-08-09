
#  DRAFT
#  DRAFT
#  DRAFT
#  DRAFT
#  DRAFT
#  DRAFT

class EMBL(object):

    _HEADERS = [
        'ACCESSION',
        'COMMENT',
        'DATE',
        'DBSOURCE',
        'DEFINITION',
        'FEATURES'
        'GENE_NAME',
        'KEYWORDS',
        'LOCUS',
        'PARENT_ACCESSION',
        'PROJECT_IDENTIFIER',
        'REFERENCE',
        'SOURCE',
    ]


    KEYS_2_SECTIONS = {
                   # identification
                   'ID': 'LOCUS',
                   'AC': 'ACCESSION',
                   # PA means PARENT ACCESSION (?) and applies to
                   # feature-level-products entries
                   'PA': 'PARENT_ACCESSION',
                   'PR': 'PROJECT_IDENTIFIER',
                   'DT': 'DATE',
                   'DE': 'DEFINITION',
                   'GN': 'GENE_NAME',  # uniprot specific
                   'KW': 'KEYWORDS',
                   # Source (taxonomy and classification)
                   'OS': 'SOURCE',
                   'OC': 'SOURCE',
                   'OG': 'SOURCE',
                   # reference keys
                   'RA': 'REFERENCE',
                   'RP': 'REFERENCE',
                   'RC': 'REFERENCE',
                   'RX': 'REFERENCE',
                   'RG': 'REFERENCE',
                   'RT': 'REFERENCE',
                   'RL': 'REFERENCE',
                   # This shuold be Reference Number. However, to split
                   # between references with _embl_yield_section I need to
                   # change section after reading one reference. So a single
                   # reference is completed when I found a new RN. The
                   # reference number information will be the reference
                   # position in the final REFERENCE list metadata
                   'RN': 'SPACER',
                   # Cross references
                   'DR': 'DBSOURCE',
                   'CC': 'COMMENT',
                   'AH': 'ASSEMBLY',
                   'AS': 'ASSEMBLY',
                   'FH': 'FEATURES',
                   'FT': 'FEATURES',
                   # sequence
                   'SQ': 'ORIGIN',
                   '  ': 'ORIGIN',
                   'CO': 'CONSTRUCTED',
                   # spacer (discarded)
                   'XX': 'SPACER'
                  }

    KEYS_TRANSLATOR = {
                   # identification
                   'ID': 'LOCUS',
                   'AC': 'ACCESSION',
                   # PA means PARENT ACCESSION (?) and applies to
                   # feature-level-products entries
                   'PA': 'PARENT_ACCESSION',
                   'PR': 'PROJECT_IDENTIFIER',
                   'DT': 'DATE',
                   'DE': 'DEFINITION',
                   'GN': 'GENE_NAME',  # uniprot specific
                   'KW': 'KEYWORDS',
                   # Source (taxonomy and classification)
                   'OS': 'ORGANISM',
                   'OC': 'taxonomy',
                   'OG': 'organelle',
                   # reference keys
                   'RA': 'AUTHORS',
                   'RP': 'REFERENCE',
                   'RC': 'REFERENCE_COMMENT',
                   'RX': 'CROSS_REFERENCE',
                   'RG': 'GROUP',
                   'RT': 'TITLE',
                   'RL': 'JOURNAL',
                   # Cross references
                   'DR': 'DBSOURCE',
                   'CC': 'COMMENT',
                   # features
                   'FH': 'FEATURES',
                   'FT': 'FEATURES',
                   'SQ': 'ORIGIN',
}

    def __init__(self, filename):
        
