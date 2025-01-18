### This script compare the relative abundance of novel vs non-novel MAGs.

# Define parameters
#assembly_list <- c("Hybrid_UHGG_Assembly_Comp50Cont5_draft2", "HybridAssembly_Comp50Cont5_draft2", "HybridAssembly_sANI98_Comp50Cont5_draft2")
assembly <- "HybridAssembly_Comp50Cont5_draft2"
#assembly <- "HybridAssembly_sANI98_Comp50Cont5_draft2"
#assembly <- "Hybrid_UHGG_Assembly_Comp50Cont5_draft2"
prefix <- assembly
#magscreen_comp_list <- c("HybridAssembly_sANI95_Comp50Cont5_vs_SPMP", "HybridAssembly_sANI98_Comp50Cont5_vs_SPMP", "HybridAssembly_sANI95_Comp50Cont5_vs_UHGG_allMAGs", "HybridAssembly_sANI98_Comp50Cont5_vs_UHGG_allMAGs")
magscreen_comp <- "HybridAssembly_sANI95_Comp50Cont5_vs_UHGG_allMAGs"
#magscreen_comp <- "HybridAssembly_sANI98_Comp50Cont5_vs_UHGG_allMAGs"
# Choose if relative abundance is calculated as mean across all samples (all) or as mean across only samples where MAG is present in (present)
abund <- "present"

# Set working directories
#counts_dir <- "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assemblies/On_raw_fastq/"
counts_dir <- "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assemblies/On_hostDecontaminated_fastqs/"
magscreen_dir <- "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/magscreen/Novel_MAGs/"
summary_dir <- paste0(counts_dir, prefix, "/summary")
#input <- list.files(summary_dir)[grepl(glob2rx('*total*strong*relAbund*revised.tsv*'), list.files(summary_dir))]
#input <- "bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund.tsv"
input <- "bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv"

print(input)

# Read in data
counts_table_raw <- read.table(paste0(counts_dir, prefix, "/summary/", input), header = TRUE, sep = "\t")
novel_MAGs <- read.table(paste0(magscreen_dir, magscreen_comp, "/new_species.txt"))$V1
novel_MAGs <- gsub(".fa", "", novel_MAGs)

# Read in total reads per sample
# setwd("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assembly comparison - mapping rates/")
# mapstats_input <- list.files()[grepl(glob2rx(paste0("*", prefix, "*strong*")), list.files())]
# print(mapstats_input)
# mapstats_raw <- read.table(mapstats_input, header = TRUE)
# mapstats <- mapstats_raw[,c("sample", "fastq_raw_reads", "kneaddata_input_reads")]

# Format counts table
colnames(counts_table_raw) <- gsub("X", "", colnames(counts_table_raw))
sample_idx <- match("ASMDAI005", colnames(counts_table_raw))
if(is.na(match("length", colnames(counts_table_raw)))) {
  col_idx <- match("Length", colnames(counts_table_raw))
  colnames(counts_table_raw)[col_idx] <- "length"
}
length_idx <- match("length", colnames(counts_table_raw))
#counts_table <- counts_table_raw[, c(1, length_idx, sample_idx:(ncol(counts_table_raw)-1))]
counts_table <- counts_table_raw[, c(sample_idx:(ncol(counts_table_raw)-2))] # sample-specific
#counts_table <- as.data.frame(counts_table_raw[, c(ncol(counts_table_raw))]) # sum across samples
if(abund == "all") {
  counts_table$mean_rel_abund <- rowMeans(counts_table)
} else if(abund == "present") {
  counts_table[counts_table==0] <- NA
  counts_table$mean_rel_abund <- rowMeans(counts_table, na.rm = TRUE)
}
relAbund_table <- cbind(counts_table_raw[,c(1,2,3)], counts_table$mean_rel_abund)
colnames(relAbund_table)[ncol(relAbund_table)] <- "mean_rel_abund"

# Partition relative abundance table into novel and non-novel partitions
if(is.na(match("original_bin", colnames(relAbund_table)))) {
  col_idx2 <- match("Original_bin_ID", colnames(relAbund_table))
  colnames(relAbund_table)[col_idx2] <- "original_bin"
}
relAbund_table_estMAGS <- relAbund_table[ ! relAbund_table$original_bin %in% novel_MAGs, ]
#relAbund_table_estMAGS_downsampled <- relAbund_table_estMAGS[sample(nrow(relAbund_table_estMAGS), length(novel_MAGs)), ]
relAbund_table_novelMAGS <- relAbund_table[ relAbund_table$original_bin %in% novel_MAGs, ]
df_non_novel <- data.frame(partition = "non_novel", mean_rel_abund = relAbund_table_estMAGS$mean_rel_abund)
df_novel <- data.frame(partition = "novel", mean_rel_abund = relAbund_table_novelMAGS$mean_rel_abund)
df_partitioned <- rbind(df_non_novel,df_novel)

# Histograms
#hist(relAbund_table_estMAGS_downsampled$mean_rel_abund, freq=TRUE, breaks=100, col=rgb(1,0,0,0.5), xlim=c(0,30))
#hist(relAbund_table_novelMAGS$mean_rel_abund, freq=TRUE, breaks=100, col=rgb(0,0,1,0.5), xlim=c(0,30), add=T)

# Boxplot
mean_values <- data.frame(group = c("novel", "non-novel"), mean = c(mean(relAbund_table_novelMAGS$mean_rel_abund), mean(relAbund_table_estMAGS$mean_rel_abund)))
mean_values$group <- as.factor(mean_values$group)
print(mean_values)
if(abund == "all") {
  ymax <- 0.002
} else if(abund == "present") {
  ymax <- 0.01
}
par(mar = c(5, 5, 3, 3))
#boxplot(df_partitioned$mean_rel_abund ~ df_partitioned$partition, horizontal = TRUE, las=2, col=rgb(0.5,0.65,0.85,0.85), pars=list(outpch = 21, outcol="grey10", outbg = "brown1", outcex=1), cex.axis = 0.5,  xlab="relative abundance", ylab="", main="")
#boxplot(df_partitioned$mean_rel_abund ~ df_partitioned$partition, horizontal = TRUE, ylim = c(0,0.002), las=2, col=rgb(0.5,0.65,0.85,0.85), pars=list(outpch = 21, outcol="grey10", outbg = "brown1", outcex=1), cex.axis = 0.5,  xlab="relative abundance", ylab="", main="MAG abundance by novelty")
boxplot(df_partitioned$mean_rel_abund ~ df_partitioned$partition, las=1, ylim = c(0,ymax), col=rgb(0.5,0.65,0.85,0.85), pars=list(outpch = 21, outcol="grey10", outbg = "brown1", outcex=1), cex.axis = 0.5,  ylab="relative abundance", xlab="", main="MAG abundance by novelty")
boxplot(df_partitioned$mean_rel_abund ~ df_partitioned$partition, las=1, col=rgb(0.5,0.65,0.85,0.85), pars=list(outpch = 21, outcol="grey10", outbg = "brown1", outcex=1), cex.axis = 0.5,  ylab="relative abundance", xlab="", main="MAG abundance by novelty")
#points(mean_values$group ~ mean_values$mean, bg = "gold", pch = 23, cex=2)
points(mean_values$mean ~ mean_values$group, bg = "gold", pch = 23, cex=2)

# Two-sample hypothesis testing. Let's use Mann-Whitney U / Wilcoxon rank-sum test, since we don't assume normal distrubution (empirically doesn't follow normal distribution)
wilcox.test(df_non_novel$mean_rel_abund, df_novel$mean_rel_abund, alternative = "two.sided")
