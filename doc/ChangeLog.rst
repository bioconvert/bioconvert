Whats' new, what has changed
================================




:Revision 0.4.4: 11 March 2020

  - BUG FIXES:

    - Fix wrong infered output names in batch mode (https://github.com/bioconvert/bioconvert/issues/252)
  - NEW converters:

    -add fasta2fasta_agp module


:Revision 0.4.3: 24 Aug 2019

 - Fix version in the setup which did not match the releas
 - Remove deprecated warning
 - Fix syntax error in embl.py module


:Revision 0.4.2: 21 Aug 2019

  - BUG Fixes:

    - Correct output extension for implicit output when input is compressed
      (https://github.com/bioconvert/bioconvert/issues/250)


:Revision 0.4.1: 12th Aug 2019

  - BUG FIXES:

    - nexus2fasta fix bug #235

:Revisin 0.4.0: (Aug 2019)

  - NEW:

    - 1-to-many and many-to-1 converters are now part of the API (Sulyvan
      Dollin)
    - Complete documentation with formats and bibliography
    - Validation of all currently implemented converters
    - New layout for the output graph
  - NEW  converters:

    - cram2fasta
    - fasta2faa
  - CHANGES:

    - cram2sam, bam2cram and related cram files have now a --reference option.
    - removed fastq2bam and fasta2bed (no clear output format)
  - BUG FXIES:

    - closed all bugs not related to new converters implementation (up to #234)

:Revision 0.3.1: (May/June 2019)

  - NEW converters:
    - bigbed2bed based on pyBigWig
    - abi2qual
  - NEW methods:
    - new method in bigwig2bedgraph based on pyBigWig
    - new method in twobit2fasta based on py2bit
  - NEW features:
    - new module io/sniffer.py
  - NEW bioconvert behaviour:
    - bioconvert --threads option added
    - bioconvert --extra-arguments option added
  - CHANGES:
    - rename the ggf3gff2 converter correctly into gff32gff2 and add test
    - fastq2fasta: (i) remove quality-file option to export qual (ii) remove
      method python_external (issues #215)
    - bam2bedgraph changes: now also reports 0-coverage bases.
    - bioconvert has now a --level option (same as --verbosity)
    - bioconvert -c/--methods is now set to -m/--methods and 
  - OTHERS:
    - update setup.py to include all scripts in ./misc package (#214)
    - More features in bioconvert_stats standalone
    - update the singularity recipes for v 0.3.0
    - update bioconda recipes (on bioconda-recipes)
    - Fix https://github.com/bioconvert/bioconvert/issues/204
    - add pytest-env plugin in the developer requirements
    - mv utils/ modules  into new directory ./io
  - BUG FIXES:
    - https://github.com/bioconvert/bioconvert/issues/226

:Revision 0.3.0: April/May 2019

    - refactoring of the core including the implementation of an implicit
      conversion. So, we can use bioconvert without specifying the conversion
      when there is no ambiguity (based on the extensions provided).
    - New working version for pypi.
    - new converters: wig2bed
    - new methods in various converters related to  phylogeny
    - Fixing the goalign and gotree scripts 
    - Fixing all tests on Travis

:Revision 0.2.0: 7 Aug 2018

    - remove pandoc from requirements. The version installed is the one from
      conda. However, conda version is a standalone only, not the code source
      so, using it in the requirements causes trouble on bioconda when building
      bioconvert. We do not use pandoc in the code for now, so let us remove it.
    - added tests for gff22gff3; moved test data sets related to GFF2 and GFF3.

:Revision 0.1.3: 2nd Aug 2018

    - add abi2fastq
    - add abi2fasta
    - add script called bioconvert_stats
    - refactorise scf2fasta and sc2fastq putting common code in utils.scf
    - update go version in requirements + update of modules using go

:Revision 0.1.1:

    - add phylip2xmfa and xmfa2phylip,
    - add bedgraph2bed
    - Fixed #132 (bedgraph2bigwig can now have chrom sizes as input)
    - add bam2wiggle, bed2wiggle, bcf2wiggle, vcf2wiggle etc

:Revision 0.1.0: 3 April 2018

    Major refactoring to allow sub commands to be used::

        bioconvert fastq2fasta test.fastq test.fasta

    instead of::

        bioconvert test.fastq test.fasta

    as well as transitive conversion: if a conversion is not implemented but
    a path exists, then conversion can be performed using -a option::

        bioconvert A2C test.A test.C -a

    - new converters: maf2sam, ods2csv, xls2csv, xlsx2csv


:Revision 0.0.12:

    - new converters: embl2fasta, embl2genbank, fasta2twobit, twobit2fasta
      fasta2fasta, sra2fastq
    - refactoring of the extensions framework to simplify the code


:Revision 0.0.11:

     - add abiliy to use go executables (add go to the requirements)
     - added converters: fasta2nexux, newick2nexus, newick2phyloxml,
       nexus2fasta, nexus2newick, nexus2phylip, nexus2phyloxml, phylip2nexus,
       phyloxml2newick, phyloxml2nexus, genbank2embl, genbank2fasta,
       stockholm2clustal and clustal2phylip

:Revision 0.0.10:

    - added samlint validator

:Revision 0.0.9:

    - added sam2paf

:Revision 0.0.8:

    - added compressor decorator
    - update bioconvert main script with several options 
    - new converters: dsrc2gz, bam2bigwig draft
    - provided squizz on bioconda and added as dependencies
    - added paflint validator

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
