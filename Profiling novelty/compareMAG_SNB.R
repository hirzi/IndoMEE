### This script compares the social niche breadths (von Meijenfeldt et al.2023; https://github.com/MGXlab/social_niche_breadth_SNB) of novel vs non-novel MAGs.

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
#colnames(counts_table)[3] <- "taxonomic lineage"

# Merge duplicate entries (add duplicates' counts)
# counts_table_dedup <- aggregate(counts_table[colnames(counts_table)[2:length(colnames(counts_table))]], by=counts_table['classification'], sum)
# counts_table_dedup <- counts_table_dedup[order(counts_table_dedup$classification),]

# Annotate (deduplicated) SNB outputwith MAG names
SNB_raw <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Social niche breadth/SNB_matrixOutput_1_.txt", header = TRUE, sep = "\t")
colnames(SNB_raw)[1] <- 'classification'
bin_class_crossref <- as.data.frame(cbind(counts_table_raw$original_bin, counts_table$classification))
colnames(bin_class_crossref) <- c("original_bin", "classification")
SNB_annot <- merge(bin_class_crossref, SNB_raw, by = 'classification')

# Load data - novel and known MAGs
magscreen.path <- "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/magscreen/Novel_MAGs/"
magscreen_comp <- "HybridAssembly_sANI95_Comp50Cont5_vs_UHGG_allMAGs"
novel.MAGs <- scan(paste0(magscreen.path, magscreen_comp, "/new_species.txt"), what = "")
novel.MAGs <- gsub(".fa", "", novel.MAGs)
all.MAGs <- scan("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/Functional enrichment/dereplicated_genomes.txt", what = "")
all.MAGs <- gsub(".fa", "", all.MAGs)
known.MAGs <- all.MAGs[which(!all.MAGs %in% novel.MAGs)]

# Partition SNB output into novel and known sets
SNB.novel.MAGs <- SNB_annot[which(SNB_annot$original_bin %in% novel.MAGs),]
SNB.known.MAGs <- SNB_annot[which(SNB_annot$original_bin %in% known.MAGs),]

# Remove duplicate MAGs present in both novel and known sets
dup_MAGs_in_novel_known <- merge(SNB.novel.MAGs, SNB.known.MAGs, by = 'classification')
SNB.novel.MAGs.dedup <- SNB.novel.MAGs[which(!SNB.novel.MAGs$classification %in% dup_MAGs_in_novel_known$classification),]
SNB.known.MAGs.dedup <- SNB.known.MAGs[which(!SNB.known.MAGs$classification %in% dup_MAGs_in_novel_known$classification),]
SNB.novel.MAGs.dedup$class <- "novel"
SNB.known.MAGs.dedup$class <- "known"

# Plot boxplot
SNB_novel_bp <- SNB.novel.MAGs.dedup[,c(5,6)]
SNB_known_bp <- SNB.known.MAGs.dedup[,c(5,6)]
SNB_bp <- rbind(SNB_novel_bp, SNB_known_bp)
mean_values <- data.frame(group = c("novel", "non-novel"), mean = c(mean(SNB_novel_bp$SNB.score, na.rm = TRUE), mean(SNB_known_bp$SNB.score, na.rm = TRUE)))
mean_values$group <- as.factor(mean_values$group)
par(mar = c(5, 5, 3, 3))
boxplot(SNB_bp$SNB.score ~ SNB_bp$class, horizontal = FALSE, las=2, col=rgb(0.5,0.65,0.85,0.85), pars=list(outpch = 21, outcol="grey10", outbg = "brown1", outcex=1), cex.axis = 0.5,  xlab="SNB score", ylab="", main="SNB score by novelty")
points(mean_values$mean ~ mean_values$group, bg = "gold", pch = 23, cex=2)

# Two-sample hypothesis testing. Let's use Mann-Whitney U / Wilcoxon rank-sum test, since we don't assume normal distrubution (empirically doesn't follow normal distribution)
wilcox.test(SNB_known_bp$SNB.score, SNB_novel_bp$SNB.score, alternative = "two.sided")
