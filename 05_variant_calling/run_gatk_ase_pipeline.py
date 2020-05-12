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

# Specify the overhang size for genome indexing by STAR.
# The overhang is maximum read length minus one.
OVERHANG      = 149

# Specify the custom index of Bowtie2-build to use
INDEX         = "Stru_genome"

# Specify the number of threads
MAX_THREADS   = 10

# Specify coverage file for read depth counts
# genomeCoverageBed requires a genome coverage file which is a tab-delimited file
# with two columns: column one: chromosome name; column two: chromosome length
# The file can be generated by running 'samtools faidx <genome.fasta>', which will
# generate an indexed fasta file with the first column containing the chromosome
# name and the second column containing chromosome lenth. The two columns can be
# cut (cut -f 1,2 <file>) into a genome index file like the one below.
COVERAGE_FILE = "data/Salmo_trutta.fSalTru1.1.dna.toplevel_genomic_index.txt"

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
INDEXED_DATA       = ["database"]

# Assign path for clean read data for aligning
RAW_DATA           = ["input/{samples}_1.fq.gz",
                      "input/{samples}_2.fq.gz"]

#def get_fq1(wildcards):
#    # code that returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
#    return sorted(glob.glob('input/' + {samples} + '*_1.fq.gz'))

#def get_fq2(wildcards):
#    # code that returns a list of fastq files for read 2 based on *wildcards.sample*, e.g.
#    return sorted(glob.glob('input/' + {samples} + '*_2.fq.gz'))

# Output aligned data and unmapped data in these directories, respectively
FIRST_PASS_ALIGN      = "results/01_first_pass/combined_Aligned.out.sam"

# Output novel splice junctions here:
SPLICE_JUNCTIONS      = "results/01_first_pass/combined_SJ.out.tab"

# Output second pass index per sample:
SECOND_PASS_INDEX     = ["results/02_second_pass_index/"]

# Output splice-aware alignments here:
SECOND_PASS_ALIGN     = "results/03_second_pass_align/{samples}Aligned.out.sam"

# Output sorted and merged BAM files here
SORTED_DATA           = "results/04_second_pass_align_sorted/{samples}_sorted_RGadded.bam"

# Output sorted and merged BAM files here
MARKED_DATA           = "results/05_second_pass_align_dup_marked/{samples}_sorted_RGadded_dedupped.bam"

# Output deduplicated metrics here:
MARKED_METRICS        = "results/06_second_pass_align_dup_metrics/{samples}_metrics.txt"

# Output split data here:
SPLIT_DATA            = "results/07_second_pass_align_dup_split/{samples}_rg_added_sorted_dedupped_split.bam"

RECALIBRATION_REPORT  = "results/08_base_recalibration/{samples}_recalibration_report.grp"

RECALIBRATED_DATA     = "results/08_base_recalibration/{samples}_rg_added_sorted_dedupped_split_recal.bam"

VARIANT_CALLS         = "results/08_variant_calling/{samples}.vcf"

FILTER_CALLS          = "results/09_filter_variant_calling/{samples}.vcf"

CREATE_VCF_LIST       = "results/10_vcf_list/combined_list.list"

COMBINED_VCF          = "results/10_vcf_list/combined.vcf.gz"

SELECT_CALLS          = "results/10_biallelic_snps/{samples}_biallelic_snps.vcf"

ASE_CALLS             = "results/11_ase_output/{samples}_ase.csv"

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
tools['gatk']    = os.path.join(dirs['src'], 'gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar')

##############################################################################
#
# Specify rules with commands to be executed
#
##############################################################################
# First rule is list the final output
rule all:
    input: expand(ASE_CALLS, samples=SAMPLES)

# Each rule will specify an intermediate step
# Align raw sequences against indexed genome
rule star_align_to_genome_first_pass:
    input:
          fq1= expand("input/{samples}_1.fq.gz", samples=SAMPLES),
          fq2= expand("input/{samples}_2.fq.gz", samples=SAMPLES)
    output: FIRST_PASS_ALIGN, SPLICE_JUNCTIONS
    threads: MAX_THREADS
    run:
        align_list=re.split('A', (''.join(output[0])))[0]
        print(input.fq1)
        input_str_fq1 = ",".join(input.fq1)
        input_str_fq2 = ",".join(input.fq2) if input.fq2 is not None else ""
        input_str =  " ".join([input_str_fq1, input_str_fq2])
	print(input_str)
	print(align_list)
        check_files_arent_empty(input)
        shell("{tools[STAR]} --readFilesCommand zcat \
                             --genomeDir {INDEXED_DATA} \
                             --readFilesIn {input_str} \
                             --outFileNamePrefix {align_list} \
                             --runThreadN {threads} \
                             --limitOutSJcollapsed 5000000 \
                             --limitIObufferSize 300000000 && \
                             [[ -s {output[0]} ]]")

# Create second index using splice junctions:
rule create_splice_aware_index:
    input: SPLICE_JUNCTIONS
    output: directory(SECOND_PASS_INDEX)
    threads: MAX_THREADS
    run:
       print(input)
       print(output)
       check_files_arent_empty(input)
       shell("{tools[STAR]} --runMode genomeGenerate \
                            --genomeDir {output} \
                            --genomeFastaFiles {GENOME} \
                            --sjdbFileChrStartEnd {input[0]} \
                            --limitSjdbInsertNsj 2000000 \
                            --sjdbOverhang 149 \
                            --runThreadN {threads}")

# Align raw sequences against indexed genome
rule star_align_to_genome_second_pass:
    input: SECOND_PASS_INDEX,
           RAW_DATA,
    output: SECOND_PASS_ALIGN
    threads: MAX_THREADS
    run:
        align_list=re.split('A', (''.join(output[0])))[0]
        print(align_list)
        check_files_arent_empty(input)
        shell("{tools[STAR]} --readFilesCommand zcat \
                             --genomeDir {input[0]} \
                             --readFilesIn {input[1]} {input[2]} \
                             --outFileNamePrefix {align_list} \
                             --sjdbOverhang 149 \
                             --runThreadN {threads} && \
                             [[ -s {output[0]} ]]")

# Align raw sequences against indexed genome
rule add_read_group_sort_compress:
    input: SECOND_PASS_ALIGN
    output: SORTED_DATA
    threads: MAX_THREADS
    run:
        sample_name=re.split('/|A', (''.join(output[0])))[2]
        check_files_arent_empty(input)
        shell("java -jar {tools[picard]} \
               AddOrReplaceReadGroups \
               I={input} \
               O={output} \
               SO=coordinate \
               RGID={sample_name} \
               RGLB=nebnext \
               RGPL=illumina \
               RGPU=novaseq6000 \
               RGSM={sample_name} && \
               [[ -s {output[0]} ]]")

# Align raw sequences against indexed genome
rule mark_duplicates:
    input: SORTED_DATA
    output: MARKED_DATA, MARKED_METRICS
    threads: MAX_THREADS
    run:
        check_files_arent_empty(input)
        shell("java -jar {tools[picard]} \
               MarkDuplicates \
               I={input} \
               O={output[0]} \
               CREATE_INDEX=true \
               VALIDATION_STRINGENCY=SILENT \
               M={output[1]} && \
               [[ -s {output[0]} ]]")

# Align raw sequences against indexed genome
rule create_genome_index:
    input: GENOME
    output: GENOME_DICT
    threads: MAX_THREADS
    run:
        check_files_arent_empty(input)
        shell("java -jar {tools[picard]} \
               CreateSequenceDictionary \
               REFERENCE={input} \
               OUTPUT=={output} && \
               [[ -s {output[0]} ]]")

# Align raw sequences against indexed genome
rule split_cigars:
    input: GENOME, MARKED_DATA
    output: SPLIT_DATA
    run:
        check_files_arent_empty(input)
        shell("java -jar {tools[gatk]} \
               SplitNCigarReads \
               -R={input[0]} \
               -I={input[1]} \
               -O={output} && \
               [[ -s {output[0]} ]]")

# Align raw sequences against indexed genome
rule variant_calls:
    input: GENOME, SPLIT_DATA
    output: VARIANT_CALLS
    run:
        check_files_arent_empty(input)
        shell("java -jar {tools[gatk]} \
               HaplotypeCaller \
               -R={input[0]} \
               -I={input[1]} \
	           --dont-use-soft-clipped-bases true \
               --standard-min-confidence-threshold-for-calling 20.0 \
               -O={output} && \
               [[ -s {output[0]} ]]")


# Align raw sequences against indexed genome
rule filter_variant_calls:
    input: GENOME, VARIANT_CALLS
    output: FILTER_CALLS
    run:
        check_files_arent_empty(input)
        shell("""java -jar {tools[gatk]} \
               VariantFiltration \
               -R={input[0]} \
               -V={input[1]} \
               -window 35 \
               -cluster 3  \
               -filter-name FS \
               -filter-expression "FS < 0.3" \
               -filter-name QD \
                -filter-expression "QD < 2.0" \
               -O={output} && \
               [[ -s {output} ]]""")

# Align raw sequences against indexed genome
rule select_variant_calls:
    input: GENOME, FILTER_CALLS
    output: SELECT_CALLS
    run:
        check_files_arent_empty(input)
        shell("""java -jar {tools[gatk]} \
            SelectVariants \
            -R={input[0]} \
            -V={input[1]} \
            -select "QD > 10.0" \
            -select 'vc.isNotFiltered()'  \
            --restrict-alleles-to BIALLELIC \
            --select-type-to-include SNP \
            -O={output} && \
            [[ -s {output} ]]""")

# Align raw sequences against indexed genome
rule ase_calls:
    input: GENOME, SPLIT_DATA, SELECT_CALLS
    output: ASE_CALLS
    run:
        check_files_arent_empty(input)
        shell("""
              java -jar {tools[gatk]} \
               ASEReadCounter \
               -R={input[0]} \
               -I={input[1]} \
               -V {input[2]} \
               --output-format CSV \
               -O={output}
               """)

# Align raw sequences against indexed genome
#rule create_vcf_list:
#    input: expand(FILTER_CALLS, samples=SAMPLES)
#    output: CREATE_VCF_LIST
#    run:
#        check_files_arent_empty(input)
#        shell("echo {input} >> {output} && \
#               [[ -s {output} ]]")


# Align raw sequences against indexed genome
#rule merge_calls:
#    input: CREATE_VCF_LIST
#    output: COMBINED_VCF
#    run:
#        check_files_arent_empty(input)
#        shell("java -jar {tools[picard]} \
#               MergeVcfs \
#               I={input} \
#               O={output} && \
#               [[ -s {output} ]]")

#java -jar GenomeAnalysisTK.jar -T VariantFiltration -R hg_19.fasta -V input.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o output.vcf 
# Align raw sequences against indexed genome
#rule base_recalibration:
#    input: GENOME, SPLIT_DATA
#    output: RECALIBRATION_REPORT, RECALIBRATED_DATA
#    run:
#        check_files_arent_empty(input)
#        shell("java -jar {tools[gatk]} \
#               PrintReads \
#               -R={input[0]} \
#               -I={input[1]} \
#               -BQSR  {output[0]} \
#               -O={output[1]} && \
#               [[ -s {output[1]} ]]")
# Convert SAM file to BAM format and sort
#rule star_align_to_genome_second_pass:
#    input: ALIGNED_DATA
#    output: SORTED_DATA
#    run:
#        check_files_arent_empty(input)
#        shell("{tools[samtools]} view -bS {input} \
#             | {tools[samtools]} sort - -f {output} && [[ -s {output} ]]")

## Add read group information for running freebayes:
#rule add_RG_information:
#    input: SORTED_DATA
#    output: RG_ANNOTATED_DATA
#    run:
#        check_files_arent_empty(input)
#        shell("java -jar {tools[picard]} AddOrReplaceReadGroups \
#                I={input} \
#                O={output} \
#                RGID={input} \
#                RGLB=lib1 \
#                RGPL=illumina \
#                RGPU=unit1 \
#                RGSM={input} \
#                && [[ -s {output} ]]")

# Calculate read depth per genomic scaffold base and coverage percentage across individual genomic
# scaffolds
#rule calculate_coverage:
#    input:  SORTED_DATA
#    output: COVERAGE_DATA
#    run:
#        check_files_arent_empty(input)
#        shell("{tools[coverage]} -ibam {input} -g {COVERAGE_FILE} > {output} && [[ -s {output} ]]")
