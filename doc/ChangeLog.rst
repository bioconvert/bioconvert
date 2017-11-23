Whats' new, what has changed
================================

:Revision 0.0.7:

    - added converters: bam2cram and cram2bam, vcf2bed
    - new class: Bioconvert that could be used for all converters !
    - new modules: core/shell, core/converter
    - add singularity to provide dot executable in the documentation and remove
      the pygraphviz dependency.

:Revision 0.0.6: added converters: bcf2vcf; vcf2bcf; bam2json; gz2bz2, bz22gz,
    gz2dsrc, .... benchmarking implemented.

:Revision 0.0.5: added bioconvert_init standalone to help developers. 
                 added gz2bz2 converter. switch default of bam2fasta with
                 sambamba

:Revision 0.0.4: update requirements and MANIFEST; added fastq2fasta, gfa2fasta

:Revision 0.0.3: benchmark in place; added fastq2fasta, scf2fastq, scf2fastq

:Revision 0.0.2: setup travis, RTD, tests; added bam2sam

:Revision 0.0.1: add bioconvert tree structure; added bam2bed, json2yaml... 
