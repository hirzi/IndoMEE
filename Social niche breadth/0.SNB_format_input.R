### This script converts metamap relative abundance tables to the input format required for Social Niche Breadth (SNB) analysis. For reference, see: https://github.com/MGXlab/social_niche_breadth_SNB

# Load libraries
library(gtools)

# Set working directories
counts_dir <- "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assemblies/On_hostDecontaminated_fastqs/"
#counts_dir <- "/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metamap/metamap_on_metaWrap_MAGs_results/"

# Set input (format) parameters
first_sample <- "ASMDAI005"

# Input assemblies
#assembly_list <- c("CoAssembly_draft1", "CoAssembly_sANI98_draft2", "HybridAssembly_Comp50Cont5_draft2", "HybridAssembly_sANI98_Comp50Cont5_draft2", "Hybrid_UHGG_Assembly_Comp50Cont5_draft2", "Hybrid_UHGG_SPMP_Assembly_Comp50Cont5_draft2", "IndAssembly_draft2", "IndAssembly_draft2_sANI98", "SPMP", "UHGG")
prefix <- "HybridAssembly_Comp50Cont5_draft2"

# List of filtering regimes
filt <- '*total*strong*relAbund*revised.tsv'

# Define input (metamap counts output) files
summary_dir <- paste0(counts_dir, prefix, "/summary")
input <- list.files(summary_dir)[grepl(glob2rx(filt), list.files(summary_dir))]
print(input)

# Read in data
counts_table_raw <- read.table(paste0(counts_dir, prefix, "/summary/", input), header = TRUE, sep = "\t")
#counts_table_raw <- read.table("/Users/hl636/Desktop/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv", header = TRUE, sep = "\t")

# Reformat counts table
sample_idx <- match(first_sample, colnames(counts_table_raw))
#counts_table <- counts_table_raw[, c(1,2,3, sample_idx:ncol(counts_table_raw))]
counts_table <- counts_table_raw[, c(3, sample_idx:ncol(counts_table_raw))]
#asd <- counts_table_raw[, c(sample_idx:ncol(counts_table_raw))]
counts_table$classification <- gsub("d__", "kingdom.", counts_table$classification)
counts_table$classification <- gsub("p__", "phylum.", counts_table$classification)
counts_table$classification <- gsub("c__", "class.", counts_table$classification)
counts_table$classification <- gsub("o__", "order.", counts_table$classification)
counts_table$classification <- gsub("f__", "family.", counts_table$classification)
counts_table$classification <- gsub("f__", "family.", counts_table$classification)
counts_table$classification <- gsub("g__", "genus.", counts_table$classification)
counts_table$classification <- gsub("s__", "species.", counts_table$classification)
counts_table$classification <- gsub(" ", "_", counts_table$classification)
counts_table$classification <- gsub("kingdom.", "super kingdom.", counts_table$classification)
counts_table <- counts_table[order(counts_table$classification),]
#colnames(counts_table)[3] <- "taxonomic lineage"

# Merge duplicate entries (add duplicates' counts)
counts_table_dedup <- aggregate(counts_table[colnames(counts_table)[2:length(colnames(counts_table))]], by=counts_table['classification'], sum)
counts_table_dedup <- counts_table_dedup[order(counts_table_dedup$classification),]

# Write out table
#write.table(counts_table, paste0(summary_dir,"/",input_pref,"_relAbund.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
write.table(counts_table_dedup, "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Social niche breadth/SNB_input_1.txt", row.names = FALSE, quote = FALSE, sep = "\t")
