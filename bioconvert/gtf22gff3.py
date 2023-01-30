"""GTF2 (General Feature Format version 2) and GFF3 (General Feature Format version 3) are both file formats used to represent genomic features and their locations in a genome. However, there are some key differences between the two formats:

Columns: GTF2 files have nine columns, while GFF3 files have nine columns plus one optional column. The optional column in GFF3 files is used to specify the sequence version.
Attributes: GTF2 files use a semicolon-separated list of key-value pairs to specify attributes, while GFF3 files use a tab-separated list of key=value pairs.
Hierarchy: GTF2 files do not have a hierarchy for features, while GFF3 files have a hierarchy where each feature is a child of a parent feature.
Sequence version: GTF2 files do not have a field for specifying the sequence version, while GFF3 files have an optional field for this purpose.
Support for multiple sequences: GTF2 files do not support multiple sequences, while GFF3 files can represent multiple sequences in a single file.
"""


def gtf2_to_gff3(gtf2_file, gff3_file):
  # Open the input GTF2 file and output GFF3 file
  with open(gtf2_file, 'r') as gtf2, open(gff3_file, 'w') as gff3:
    # Iterate over each line in the GTF2 file
    for line in gtf2:
      # Skip comments in the GTF2 file
      if line.startswith('#'):
        continue

      # Split the line by tabs
      fields = line.strip().split('\t')

      # Extract the relevant fields from the GTF2 file
      seqid = fields[0]
      source = fields[1]
      feature = fields[2]
      start = fields[3]
      end = fields[4]
      score = fields[5]
      strand = fields[6]
      phase = fields[7]
      attributes = fields[8]
      
      # Convert the GTF2 attributes field to a dictionary
      attribute_dict = {}
      for attribute in attributes.split(';'):
        attribute = attribute.strip()
        if not attribute:
          continue
        key, value = attribute.split(' ', 1)
        attribute_dict[key] = value.strip('"')
      
      # Extract the ID and Name attributes from the dictionary
      # If the ID attribute is not present, use the Name attribute
      # If the Name attribute is not present, use the ID attribute
      # If neither the ID nor Name attribute is present, use '.'
      if 'ID' in attribute_dict:
        id_attribute = attribute_dict['ID']
      elif 'Name' in attribute_dict:
        id_attribute = attribute_dict['Name']
      else:
        id_attribute = '.'
      
      if 'Name' in attribute_dict:
        name_attribute = attribute_dict['Name']
      elif 'ID' in attribute_dict:
        name_attribute = attribute_dict['ID']
      else:
        name_attribute = '.'
      
      # Write the GFF3 line to the output file
      gff3.write(f"{seqid}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{phase}\tID={id_attribute};Name={name_attribute}\n")






"""



# seqid  source  feature start   end     score   strand  phase   attributes
chr1    source1 gene    11874   14409   .       +       .       gene_id "ENSG000001"; gene_name "DDX11L1";
chr1    source1 transcript      11874   14409   .       +       .       gene_id "ENSG000001"; gene_name "DDX11L1"; transcript_id "ENST00000456328";
chr1    source1 exon    11874   12227   .       +       .       gene_id "ENSG000001"; gene_name "DDX11L1"; transcript_id "ENST00000456328"; exon_number "1";
chr1    source1 exon    12613   12721   .       +       .       gene_id "ENSG000001"; gene_name "DDX11L1"; transcript_id "ENST00000456328"; exon_number "2";
chr1    source1 exon    13221   14409   .       +       .       gene_id "ENSG000001"; gene_name "DDX11L1"; transcript_id "ENST00000456328"; exon_number "3";
"""
