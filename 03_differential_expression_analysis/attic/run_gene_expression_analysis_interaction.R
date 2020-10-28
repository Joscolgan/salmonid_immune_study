#!/usr/bin/env Rscript

libraries <- c("tximport",
"DESeq2")
for (lib in libraries) {
    if (require(package = lib, character.only = TRUE)) {
        print("Successful")
    } else {
        print("Installing")
        source("https://bioconductor.org/biocLite.R")
        library(lib, character.only = TRUE )
    }
}
dir.create("results")

## Take input from the command line:
args <- commandArgs(TRUE)
input_files <- args[1] # Assign first argument as input
input_samples <- args[2]
output <- args[3]						# Assign second argument as input
message(paste("will output to", output))

## Set paths to folders containing output files from kallisto
paths <- list()
paths$sample_batch_info <- input_files
paths$salmon_output   <- input_samples

# Set relative paths:
paths$salmon_files_relative <- grep(x = list.files(paths$salmon_output,
                                                     recursive = TRUE),
                                      pattern = ".sf",
                                      value   = TRUE)
paths$salmon_files <- file.path(paths$salmon_output,
                                  paths$salmon_files_relative)

# Automatically extract file names from kallisto output for colnames
names(paths$salmon_files) <- gsub(paths$salmon_files_relative,
                                    pattern = "/.*",
                                    replacement = "")

#for (filenumber in 1:length(paths$salmon_files)) {
#  current_name <- names(paths$salmon_files)[filenumber]
#  current_file <- paths$salmon_files[filenumber]
#  if (FALSE == grepl(pattern = current_name, x = current_file)) {
#    kill("we have a problem - names and filenames dont match up")
#  }
#}

## Extract sample names and put into df for tximport
samples     <- data.frame(treatment = names(paths$salmon_files))

## Read in sample information for room and batch:
samples_information <- read.table(file = paths$sample_batch_info,
                                  header = FALSE,
                                  col.names = c("sample_name",
                                                "sex",
                                                "batch",
"family",
"population",
"group"),
                                  row.names = 1)

## Ensure family is a factor:
samples_information$family <- as.factor(samples_information$family)


## Make sure 'CON' is the reference:
samples_information$sex <- relevel(x = samples_information$sex,
                                         ref = "F")

# sanity check
#if (FALSE == all(samples$treatment %in% rownames(samples_information))) {
#    kill("we have a problem - sample names from files and
#         from batch descriptor file  dont match up")
#}
samples_information

# Read in file corresponding to transcripts to gene ids
transcript_to_gene <- read.table(file = "./data/combined_transcript_to_gene.txt",
                                 col.names = c("transcript",
                                               "locus"))
# Use tximport on kallisto files
## Counts are estimated counts from kallisto:
txi_counts <- tximport(paths$salmon_files,
                       type    = "salmon",
                       tx2gene = transcript_to_gene,
                       countsFromAbundance = "no")

# Save object for use in other scripts
dir.create(path = output)
save(txi_counts, file = paste(output,
                              "txi_count_estimates.Rdata", sep =""))

## Ensure column and row names are the same:
print(colnames(txi_counts$counts))
print(rownames(samples_information))
rownames(samples_information) <- colnames(txi_counts$counts)

## Construct DESeq data set from txi object and sample information:
deseq_txi <- DESeqDataSetFromTximport(txi     = txi_counts,
                                      colData = samples_information,
                                      design  = ~population + sex + population:sex)

## Filter out lowly expressed genes:
keep <- rowSums(counts(deseq_txi)) >= 10
deseq_txi <- deseq_txi[keep, ]
print(nrow(deseq_txi))

save(deseq_txi, file = paste(output,
                              "deseq_txi_object.Rdata", sep =""))


# Run DESeq2 (Wald test) for differential expression analysis:
## Which genes are differentially expressed within each treatment for
## each caste in comparison to control?
## Perform a pairwise comparison between treatments:
deseq_object  <- DESeq(deseq_txi,
                       test = "Wald",
                       betaPrior = FALSE)

deseq_results <- results(deseq_object,
                         contrast = c("sex",
                                      "F",
                                      "M"))
print(nrow(deseq_results))


## Convert into a dataframe:
deseq_results_df <- as.data.frame(deseq_results)
print(nrow(deseq_results_df))
# Results table ordered by smallest p value:
deseq_results_df_reordered <- deseq_results_df[order(deseq_results_df$pvalue), ]

# Subset significant genes:
deseq_results_df_reordered_sig <- subset(deseq_results_df_reordered,
                                                         padj <= 0.05)
print(nrow(deseq_results_df_reordered_sig))

# Save to output:
write.table(x = deseq_results_df_reordered,
            file = paste(output, "/differentially_expressed_genes.txt", sep = ""),
            sep = "\t",
            quote = FALSE)

## Exploring the interaction between population and sex:
# Subset significant genes:
## Compare samples based on interaction
deseq_results <- results(deseq_object,
                         contrast = c("population",
                                      "Err",
                                      "Bun"))
print(nrow(deseq_results))
## Convert into a dataframe:
deseq_results_df <- as.data.frame(deseq_results)
print(nrow(deseq_results_df))
# Results table ordered by smallest p value:
deseq_results_df_reordered <- deseq_results_df[order(deseq_results_df$pvalue), ]
# Subset significant genes:
deseq_results_df_reordered_sig <- subset(deseq_results_df_reordered,
                                                         padj <= 0.05)
print(nrow(deseq_results_df_reordered_sig))
print(c(paste("Total number of genes:",
              nrow(deseq_results_df_reordered_sig))))

## Define interaction:
interaction <- "populationErr.sexM"
## Create results table to explore interaction term:
diff_across_populations <- results(deseq_object, name = interaction)
## Convert results table into dataframe:
diff_across_populations_df <- as.data.frame(diff_across_populations)
## Reorder by significance and count number of significantly
## differentially expressed genes:
diff_across_populations_df_sorted <- diff_across_populations_df[order(diff_across_populations_df$padj), ]
## Count the number of significant genes:
diff_across_populations_df_count <- nrow(subset(diff_across_populations_df_sorted,
                                                   padj < 0.05))
## Print to console:
print(c(paste("Total number of genes:",
              diff_across_populations_df_count)))
diff_across_populations_df_genes <- subset(diff_across_populations_df_sorted,
                                           padj < 0.05)
## Overlap with genes differentially expressed between populations?
interaction_sig_gene_list <- row.names(diff_across_populations_df_genes)
population_sig_gene_list <- row.names(deseq_results_df_reordered_sig)
table(table(sort(c(interaction_sig_gene_list,
        population_sig_gene_list))))

# Save to output:
write.table(x = diff_across_populations_df_genes,
            file = paste(output, "/differentially_expressed_genes_interaction.txt",
                         sep = ""),
            sep = "\t",
            quote = FALSE)
