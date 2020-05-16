#!/usr/bin/env python3
##############################################################################
##############################################################################
# Author: Joe Colgan                   Program: run_hisat2_to_htseq.py
#
# Date: 01/10/2019
#
# Purpose:
# The purpose of this script is to take two FASTQ files (i.e. pairs) per sample
# and aligns against an index-genome. This script performs two pass alignment
# as implemented in the short-read RNA-Seq aligner STAR. The script marks
# duplicates and sorts the output bam produced by STAR using Picard.
##############################################################################
# Import modules
import os.path
import re
import glob

from helper_functions import *

# To achieve final output, the Snakefile contains custom-defined rules (see below) that
# outline commands to execute sequentially to take custom-defined input(s) and generate
# final output (as defined in rule all).

# For this specific script, rules are defined as follows:
# rule all:
#   - Defines the expected final output of the Snakefile.
# rule index_genome:
#    - A single fasta file containing genome of interest is provided as input.
#      Input fasta file is indexed by bowtie2.
# rule align_to_genome:
#   - For each sample, sequences are aligned in pairs to the bowtie2-indexed genome.
#    The rule checks if sequence headers for each input pair contain the same header information.
#    If they differ, an error is raised.
# rule sort_sam_to_bam:
#    - For each sample, aligned SAM file is converted to BAM file and sorted.
# rule calculate_coverage:
#    - For each sample, calculates read depth per base for each genomic scaffold using an input file
#      containing genomic co-ordinates for each genomic scaffold.
#    - In addition, summarises percentage of read depth across entire genome.
# rule parse_low_coverage_regions:
#    - For each input file, parses rows where the second column contains a value less than 5.
#    - Subsequently, parses rows containing 'genome' (which is summary information) and outputs.
# rule reformat_plot_data:
#    - Using a custom R script, reformat the parsed data into a format that can be plotted.
# rule combine_plot_input:
#    - Combine reformatted data for each sample.
# rule plot_stack_charts:
#    - Using the combined reformatted data for each sample, plot a stacked bar plot.

##############################################################################
# Sample information
##############################################################################
# Brown trout (Salmo trutta) males and females sampled summer 2018
# Originating from two populations that vary in their propensity to migrate to
# sea (anadromy strategy). Erriff ("Err", migratory) and Bunaveela ("Bun", resident)
# Each individual and site were assigned unique identifiers
# Example: 'D10_2018_Stru_Err_2_F_M8_brain'
# Explanation:{extraction_id}{year_collected}_{species}_{population}_{tank}_{sex}_{tube_number}_{tissue_type}
# species: Stru = Salmo trutta
# sex: M = Male; F = Female

##############################################################################
# Prior to use
##############################################################################
# To run alignment_to_coverage.py:
#  1. Download and install the following software:
#   HISAT2
#   HTSeq
#
#  2. Ensure helper_functions.py is within the same directory of the Snakefile.
# For this particular snakefile, ensure the location of the bash script(genomecov_output_parser.sh)
#  and Rscripts (reformat_for_plot_test.R and plot_stacks.R) are defined.
#
#  3. Assign global variables for use within specific rules
#     Please see section below for further information on variable to be assigned
#
#  4. Assignment of wildcards to be used within the rules
#
#  5. Define all the input and outputs of each rule with respect to the above defined wildcards
#
#  6. Each rule willl take assigned input, global variables and generate defined outputs.
#
#  7. Input data should be formatted in the context of defined wildcards: {sample}_{pair}.fq.gz
#      For example: D10_2018_Stru_Err_2_F_M8_brain_1.fq.gz
#
# 8. Make a text.file called 'sample_list.txt' and put in same directory as Snakfile.
#     Populate 'sample_list.txt' with names of samples to be analysed.
#       For example:
#       D10_2018_Stru_Err_2_F_M8_brain
#       F12_2018_Stru_Err_2_F_M8_gonad
#       H4_2018_Stru_Bun_2_F_O22_liver

##############################################################################
# Assign global variables for use in rules (see below)
##############################################################################
# Specify the path to reference annotation file
GENOME_GTF    = "data/Salmo_trutta.fSalTru1.1.99.gtf"

# Specify the number of threads
MAX_THREADS   = 10

##############################################################################
# Assignment of wildcards to be used within rules
##############################################################################
# Open file and read in contents - one sample per line
with open('samples_list.txt') as samples:
    content = samples.readlines()
    SAMPLES = [samples.rstrip('\n') for samples in content]
    print(SAMPLES)

##############################################################################
# Specify all input/output files in terms of sample wildcards
##############################################################################
# Assign path for hisat2 splice aware bam files
HISAT2_BAM         = ["input/{samples}.bam"]

# Assign path for exon-level counts extracted from alignment files:
HTSEQ_EXON_COUNTS  = "results/02_exon_level_counts/{samples}_exon_level_counts.txt"
##############################################################################
# Define binaries in context of path relative to Snakefile
##############################################################################
# binaries
# Align lines of code using 'Assign Align' using cmd+shift+p and selecting 'align'
# Create dictionaries for directories and  tools'
dirs  = {}
dirs['project'] = os.path.abspath('/workspace/home/jcolgan/')
dirs['src']     = os.path.join(dirs['project'], '.local/bin/')

# Create an empty dictionary called 'tools'
tools = {}
tools['HTSeq-count']        = os.path.join(dirs['src'], 'htseq-count')

##############################################################################
#
# Specify rules with commands to be executed
#
##############################################################################
# First rule is list the final output
rule all:
    input: expand(HTSEQ_EXON_COUNTS, samples=SAMPLES)

# Each rule will specify an intermediate step
# Align raw sequences against indexed genome
rule extract_exon_level_counts:
    input: HISAT2_BAM, GENOME_GTF
    output: HTSEQ_EXON_COUNTS
    run:
        shell("{tools[HTSeq-count]} \
                --format=bam \
                --format=bam \
                --order=name  \
                --stranded=no \
                --a=10 \
                --type=exon \
                --mode=union \
                --nonunique=none {input} \
                {input[0]} {input[1]} > {output} && \
                [[ -s {output[0]} ]]")
