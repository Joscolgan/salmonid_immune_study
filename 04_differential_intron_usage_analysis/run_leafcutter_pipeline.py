#!/usr/bin/env python3
##############################################################################
##############################################################################
# Author: Joe Colgan                   Program: run_gatk_ase_pipeline.py
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
#   STAR
#   Picard
#   GATK
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
GENOME        = "data/Salmo_trutta.fSalTru1.1.dna.toplevel.fa"

GENOME_DICT   = "data/Salmo_trutta.fSalTru1.1.dna.toplevel.dict"

# Specify the path to reference annotation file
GENOME_GTF    = "data/Salmo_trutta.fSalTru1.1.99.gtf.gz"

EXON_FILE     = "data/Salmo_trutta.fSalTru1.1.99_exons.txt.gz"

# Specify the overhang size for genome indexing by STAR.
# The overhang is maximum read length minus one.
OVERHANG      = 149

# Specify the custom index of Bowtie2-build to use
INDEX         = "Stru_genome"

# Specify the number of threads
MAX_THREADS   = 40

# Specify coverage file for read depth counts
# genomeCoverageBed requires a genome coverage file which is a tab-delimited file
# with two columns: column one: chromosome name; column two: chromosome length
# The file can be generated by running 'samtools faidx <genome.fasta>', which will
# generate an indexed fasta file with the first column containing the chromosome
# name and the second column containing chromosome lenth. The two columns can be
# cut (cut -f 1,2 <file>) into a genome index file like the one below.
COVERAGE_FILE = "./data/GCF_901001165.1_fSalTru1.1_genomic_index.txt"

##############################################################################
# Assignment of wildcards to be used within rules
##############################################################################
# Open file and read in contents - one sample per line
with open('samples_list_erriff.txt') as samples:
    content = samples.readlines()
    SAMPLES = [samples.rstrip('\n') for samples in content]
    print(SAMPLES)

# Assign pair information
PAIR          = ["_1", "_2"]

##############################################################################
# Specify all input/output files in terms of sample wildcards
##############################################################################
# Assign output data for indexed genome data
INDEXED_DATA       = ["database_star"]

# Assign path for clean read data for aligning
RAW_DATA           = ["input/{samples}_1.fq.gz",
                      "input/{samples}_2.fq.gz"]

SAMPLE_INFORMATION = "sample_information_combined_sorted_erriff.txt"

#def get_fq1(wildcards):
#    # code that returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
#    return sorted(glob.glob('input/' + {samples} + '*_1.fq.gz'))

#def get_fq2(wildcards):
#    # code that returns a list of fastq files for read 2 based on *wildcards.sample*, e.g.
#    return sorted(glob.glob('input/' + {samples} + '*_2.fq.gz'))

# Output aligned data and unmapped data in these directories, respectively
FIRST_PASS_ALIGN      = "results/01_first_pass/{samples}_Aligned.out.bam"

SPLICE_JUNCTIONS      = "results/02_splice_junctions/{samples}.junc"

SPLICE_LIST           = "results/02_splice_junctions/combined_juncfiles_erriff.txt"

INTRON_CLUSTERS       = "results/03_intron_clusters/combined_perind_numers.counts.gz"

DIE_OUTPUT            = "results/04a_differential_intron_erriff/results_cluster_significance.txt"

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
tools['STAR']        = os.path.join(dirs['src'], 'STAR/bin/Linux_x86_64/STAR')
tools['picard']      = os.path.join(dirs['src'], 'picard/build/libs/picard.jar')
tools['gatk']        = os.path.join(dirs['src'], 'gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar')
tools['bamsh']       = os.path.join(dirs['src'], 'leafcutter/scripts/bam2junc.sh')
tools['leafcutter']  = os.path.join(dirs['src'], 'leafcutter/clustering/leafcutter_cluster.py')
tools['leafcutter_ds'] = os.path.join(dirs['src'], 'leafcutter/scripts/leafcutter_ds.R')
tools['gtf_to_exons'] = os.path.join(dirs['src'], 'leafcutter/scripts/gtf_to_exons.R')
tools['Rscript']      = os.path.join(dirs['project'], 'R/bin/bin/Rscript')

##############################################################################
#
# Specify rules with commands to be executed
#
##############################################################################
# First rule is list the final output
rule all:
    input: expand(DIE_OUTPUT)

# Align raw sequences against indexed genome
rule create_star_index:
    input: GENOME, GENOME_GTF
    output: directory(INDEXED_DATA)
    threads: MAX_THREADS
    run:
        align_list=re.split('A', (''.join(output[0])))[0]
        print(align_list)
        check_files_arent_empty(input)
        shell("{tools[STAR]} --readFilesCommand zcat \
              --runMode genomeGenerate \
              --genomeDir {output} \
              --genomeFastaFiles {input[0]} \
              --sjdbGTFfile {input[1]} \
              --sjdbOverhang 149 \
              --runThreadN {threads}")

# Each rule will specify an intermediate step
# Align raw sequences against indexed genome
rule star_align_to_genome_first_pass:
    input:INDEXED_DATA,
          RAW_DATA
    output: FIRST_PASS_ALIGN
    threads: MAX_THREADS
    run:
        align_list=re.split('A', (''.join(output[0])))[0]
        #print(input.fq1)
#        input_str_fq1 = ",".join(input.fq1)
#        input_str_fq2 = ",".join(input.fq2) if input.fq2 is not None else ""
#        input_str =  " ".join([input_str_fq1, input_str_fq2])
	#print(input_str)
	#print(align_list)
        check_files_arent_empty(input)
        shell("{tools[STAR]} --readFilesCommand zcat \
                             --genomeDir {input[0]} \
                             --readFilesIn {input[1]} {input[2]} \
                             --twopassMode Basic \
                             --outSAMstrandField intronMotif \
                             --outSAMtype BAM Unsorted \
                             --outFileNamePrefix {align_list} \
                             --runThreadN 40 \
                             --limitOutSJcollapsed 5000000 \
                             --limitIObufferSize 300000000 && \
                             [[ -s {output[0]} ]]")

## Generate a list of splice juctions:
rule extract_splice_junctions:
    input: FIRST_PASS_ALIGN
    output: SPLICE_JUNCTIONS
    threads: MAX_THREADS
    run:
        check_files_arent_empty(input)
        shell("sh {tools[bamsh]} {input} {output} && \
              [[ -s {output} ]]")

## Generate a list of splice juctions:
rule create_bam_list:
    input: expand(SPLICE_JUNCTIONS, samples=SAMPLES)
    output: SPLICE_LIST
    threads: MAX_THREADS
    run:
        check_files_arent_empty(input)
        shell("""
              echo {input} | tr ' ' '\n' > {output} && \
              [[ -s {output} ]]""")

## Cluster by introns:
rule cluster_introns:
    input:  SPLICE_LIST
    output: INTRON_CLUSTERS
    threads: MAX_THREADS
    run:
        align_list=re.split('p', (''.join(output[0])))[0]
        align_list=align_list[:-1]
        check_files_arent_empty(input)
        shell("""
              python2.7 {tools[leafcutter]} -j {input} \
              -m 50 \
              -l 500000 \
              -o {align_list} """)

## Run differential exon usage analysis:
rule run_exon_analysis:
    input:  EXON_FILE, INTRON_CLUSTERS, SAMPLE_INFORMATION
    output:  DIE_OUTPUT
    threads: MAX_THREADS
    run:
        die_output=re.split('c', (''.join(output[0])))[0]
        die_output=die_output[:-1]
        check_files_arent_empty(input)
        shell("""
              {tools[Rscript]} {tools[leafcutter_ds]} \
              --num_threads {threads} \
              --exon_file={input[0]} \
              -i 5 -g 5 \
              -p {threads} \
              {input[1]} \
              {input[2]}  \
              -o {die_output} && \
              [[ -s {output} ]]""")
