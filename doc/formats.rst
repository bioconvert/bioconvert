

BAM format
================

:reference: http://samtools.github.io/hts-specs/SAMv1.pdf



BED format
====================

:reference: http://genome.ucsc.edu/FAQ/FAQformat.html#format1

BED file must has at least 3 columns (chrom, start, end).

Wiggle Track format (WIG)
============================

:reference: http://genome.ucsc.edu/goldenPath/help/wiggle.html

The bigWig format is used for graphing track needs. The wiggle (WIG) format is
an older format for display of dense, continuous data such as GC percent. 
Wiggle data elements must be equally sized. 

Similar format such as the bedGraph format is also an older format used to display sparse data
or data that contains elements of varying size.

For speed and efficiency, wiggle data is compressed with a minor loss of precision when
data is exported from a wiggle track.
