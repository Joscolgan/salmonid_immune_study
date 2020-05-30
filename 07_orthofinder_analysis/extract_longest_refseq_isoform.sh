#!/usr/bin/env bash
####################################################################
##
## Author: Joe Colgan                  Date: 12-04-2020
##
## Program: extract_longest_refseq_isoform.sh
##
## Purpose:
## This script takes two FASTA files as input:
## 1) An uncompressed RefSeq protein FASTA file.
## 2) An uncompressed RefSeq translated cds FASTA file.
## The script takes these input as arguments from the command line.
## The script outputs a FASTA file containing a single longest
## protein isoform per gene.
####################################################################

## Take arguments from the command line:
protein_file=$1
translated_cds_file=$2

## Checked arguments exist:
if [[ $# -eq 0 ]] ; then
    echo 'No arguments supplied'
    echo 'extract_longest_refseq_isoform.sh <protein_file> <translated_cds>'
    exit 1
fi

## Create variable for sample name:
sample="$(echo "$protein_file" | cut -d '_' -f 1-3)"

## Create tmp_directory:
tmp_folder="$(echo "$sample"_tmp)"
echo "$tmp_folder"
mkdir "$tmp_folder"

## Print to console:
echo "Running script"

##
zcat "$translated_cds_file" | cut -d '[' -f 1-3 | sed 's/ /_/g' \
> $tmp_folder/"$sample"_cds_edited.fa 

## Extract gene information and create a unique list of gene ids:
grep '>' $tmp_folder/"$sample"_cds_edited.fa | \
cut -d '[' -f 3 | sed 's/db_xref=GeneID://g' | \
sed 's/]_//g' | sort | uniq > $tmp_folder/unique_gene_list.txt

## Calculate the length for each predicted protein isoform:
seqtk comp $tmp_folder/"$sample"_cds_edited.fa > $tmp_folder/length_per_translated_cds.txt

## Match lines containing gene name, sort by descending length and extract the longest isoform per gene:
while read line;
do
echo "$line";
grep -w "$line" $tmp_folder/length_per_translated_cds.txt | \
sort -k2,2nr | head -n 1 | cut -f 1  >> $tmp_folder/longest_transcript_per_gene.txt;
done < $tmp_folder/unique_gene_list.txt

## Extract protein ID:
cut -d '_' -f 4,5 $tmp_folder/longest_transcript_per_gene.txt \
> $tmp_folder/longest_transcript_per_gene_protein_ids.txt

## Extract FASTA headers of corresponding protein IDs in the RefSeq protein sequence FASTA
## The protein sequences and translated cds should be the same but want to ensure:
zgrep -f $tmp_folder/longest_transcript_per_gene_protein_ids.txt "$protein_file" | \
sed 's/>//g' > $tmp_folder/"$sample"_longest_protein_ids.txt

## The last step is to subset longest protein isoform per gene:
seqtk subseq "$protein_file" $tmp_folder/"$sample"_longest_protein_ids.txt -l 80 \
> "$sample"_primary_transcripts.fa

## Print to console:
echo "Complete!"
