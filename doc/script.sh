#!/bin/bash
FILES=*fastq
CONVERSION=fastq2fasta
for f in $FILES
do
  echo "Processing $f file..."
  bioconvert $CONVERSION $f  --force
done
