#!/bin/bash
for file in *_dedup.fastq
do
  echo "Running FastQC on $file..."
  fastqc -o ../FastQC -f fastq -t 10 "$file"
done