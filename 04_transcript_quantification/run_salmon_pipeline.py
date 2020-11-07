#!/usr/bin/env python3
##############################################################################
##############################################################################
# Author: Joe Colgan                   Program: run_salmon_pipeline.py
#
# Date: 04/10/2019
#
# Purpose:
# The purpose of this script is to take two FASTQ files (i.e. pairs) per sample
# and aligns against an index-genome. This script performs pseudoalignment using
# Salmon. The script outputs tab-dimilited text files for each sample.
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
# rule create_salmon_index
#   - Build index for user-provided transcriptome/cdna FASTA file
# rule generate_salmon_quants
#   - Defines input and output files for pseudoaligner to take input FASTQ files
#     and output transcript abundances per each sample.

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
#   Salmon
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
# Specify the reference genome file to be indexed
RNA_FASTA        = "data/gentrome.fa.gz"

# Specifi path to decoy transcriptome:
DECOYS           = "./decoys.txt"

# Specify the number of threads
MAX_THREADS   = 15

##############################################################################
# Assignment of wildcards to be used within rules
##############################################################################
# Open file and read in contents - one sample per line
with open('samples_list.txt') as samples:
    content = samples.readlines()
    SAMPLES = [samples.rstrip('\n') for samples in content]
    print(SAMPLES)

# Assign pair information
PAIR          = ["_1", "_2"]

##############################################################################
# Specify all input/output files in terms of sample wildcards
##############################################################################
# Assign output data for indexed genome data
SALMON_DATABASE       = ["database"]

# Assign path for clean read data for aligning
RAW_DATA           = ["input/{samples}_1.fq.gz",
                      "input/{samples}_2.fq.gz"]

# Output aligned data and unmapped data in these directories, respectively
SALMON_QUANTS      = "results/01_salmon_quants/{samples}/"

##############################################################################
# Define binaries in context of path relative to Snakefile
##############################################################################
# binaries
# Align lines of code using 'Assign Align' using cmd+shift+p and selecting 'align'
# Create dictionaries for directories and  tools'
dirs  = {}
dirs['project'] = os.path.abspath('/workspace/home/jcolgan/')
dirs['src']     = os.path.join(dirs['project'], 'src')

# Create an empty dictionary called 'tools'
tools = {}
tools['salmon']        = os.path.join(dirs['src'], 'salmon-latest_linux_x86_64/bin/salmon')

##############################################################################
#
# Specify rules with commands to be executed
#
##############################################################################
# First rule is list the final output
rule all:
    input: expand(SALMON_QUANTS, samples=SAMPLES)

# Each rule will specify an intermediate step
# Align raw sequences against indexed genome
rule create_salmon_index:
    input: RNA_FASTA, DECOYS
    output: directory(SALMON_DATABASE)
    threads: MAX_THREADS
    run:
        check_files_arent_empty(input)
        shell("{tools[salmon]} index \
                             -t {input[0]} \
                             -d {input[1]} \
                             -p 15 \
                             -i {output} && \
                             [[ -s {output[0]} ]]")

# Create second index using splice junctions:
rule generate_salmon_quants:
    input: SALMON_DATABASE, RAW_DATA
    output: directory(SALMON_QUANTS)
    threads: MAX_THREADS
    run:
       print(input)
       print(output)
       check_files_arent_empty(input)
       shell("{tools[salmon]} quant \
                            -i {input[0]} \
                            -l A \
                            -1 {input[1]} \
                            -2 {input[2]} \
                            -p {threads} \
                            --validateMappings \
                            -o {output} && \
                            [[ -s {output[0]} ]]")
