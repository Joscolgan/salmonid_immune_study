#!/bin/sh

## Take inputs from the command line:
input_fasta=$1
input_gtf=$2
overhang=$3

## Check arguments are provided:
if [ $# -eq 0 ]
  then
    echo "No arguments supplied. Usage: ./make_star_index.sh input.fasta input.gtf overhang"
    echo "For overhang, a value of the maximum read length minus 1 should be provided."

fi

## Create STAR indices
STAR --runThreadN 20 \
     --runMode genomeGenerate \
     --genomeDir ./database/ \
     --genomeFastaFiles "$input_fasta" \
     --sjdbGTFfile "$input_gtf" \
     --sjdbOverhang "$overhang"
