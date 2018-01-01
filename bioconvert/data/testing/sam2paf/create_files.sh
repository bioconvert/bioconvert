


#wget Hm2_GTGAAA_L005_R1_001.fastq.gz select only 100 first alignments

make clean
# this file is also in ../../
wget "https://www.ebi.ac.uk/ena/data/view/K01711&display=fasta" -O measles.fa
unpigz  -c ../../measles_R1.fastq.gz | head -n 400 > measles_R1.fastq
pigz measles_R1.fastq


#minimap2 measles.fa measles_R1.fastq.gz > approx-mapping.paf

#
minimap2 -c measles.fa measles_R1.fastq.gz > test_sam2paf_v1.paf

# or to output alignments in the SAM format:
minimap2 -a measles.fa measles_R1.fastq.gz > test_sam2paf_v1.sam

# bioconvert alignment.sam test.paf
