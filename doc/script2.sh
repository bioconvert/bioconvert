#!/bin/bash
FILES=*fastq
CONVERSION=fastq2fasta
for f in $FILES
do
  echo "Processing $f file..."
  sbatch -c 1 bioconvert $CONVERSION $f  --force
done
