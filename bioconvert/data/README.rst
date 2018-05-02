CONVENTION
---------------
test_<convertNAME>_vXX.<input_extension> 
test_vcf2bcf_v1.vcf


bcf2vcf and vcf2bcf
--------------------

test_vcf2bcf_v1.vcf  origin unknown
test_vcf2bcf_v1.bcf  created using bcftools view -O b test_vcf2bcf_v1.vcf


bam2wiggle
------------

bam input is test_measles.sorted.bam

bcf2wiggle
--------------

bcf input is test_bcf2vcf_v1.bcf


json2yaml and yaml2json
-------------------------

- test files: test_v1.json and test_v1.yaml 
- description: simple dictionary-like data sets with comments in the YAML file


gfa2fasta
------------

- test_v1.gfa reference: http://seqanswers.com/forums/showthread.php?t=64862


asqg
--------
- test_example1.asqg: from https://github.com/jts/sga/wiki/ASQG-Format documentation

- test_asqg_to_gfa.asqg: from https://github.com/sjackman/assembly-graph/blob/master/sample.gfa




PAF, GFA
-----------------

For a large example, we can use the example from
https://github.com/lh3/miniasm/blob/master/misc/demo-ecoli-pacbio.sh

Note that we must use minimap (not minimap2 for which utg.gfa is empty)::

    wget -O- http://www.cbcb.umd.edu/software/PBcR/data/selfSampleData.tar.gz | tar zxf -
    ln -s selfSampleData/pacbio_filtered.fastq reads.fq
    minimap -Sw5  -L100 -m0 -t8 reads.fq reads.fq| gzip -1 > reads.paf.gz 
    miniasm -f reads.fq reads.paf.gz > utg.gfa




