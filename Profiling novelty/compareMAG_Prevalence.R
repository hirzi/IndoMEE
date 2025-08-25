### This script compares the relative abundance and prevalence of novel vs known MAGs

# Load libraries
library(gtools)
library(matrixStats)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(stringr)

# Set working directories
counts_dir <- "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assemblies/On_hostDecontaminated_fastqs/"
magscreen_dir <- "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/magscreen/Novel_MAGs/"

# Define parameters
#assembly_list <- c("CoAssembly_draft1", "CoAssembly_sANI98_draft2", "HybridAssembly_Comp50Cont5_draft2", "HybridAssembly_sANI98_Comp50Cont5_draft2", "Hybrid_UHGG_Assembly_Comp50Cont5_draft2", "Hybrid_UHGG_SPMP_Assembly_Comp50Cont5_draft2", "IndAssembly_draft2", "IndAssembly_draft2_sANI98", "SPMP", "UHGG")
prefix <- "HybridAssembly_Comp50Cont5_draft2"
#magscreen_comp_list <- c("HybridAssembly_sANI95_Comp50Cont5_vs_SPMP", "HybridAssembly_sANI98_Comp50Cont5_vs_SPMP", "HybridAssembly_sANI95_Comp50Cont5_vs_UHGG_allMAGs", "HybridAssembly_sANI98_Comp50Cont5_vs_UHGG_allMAGs")
magscreen_comp <- "HybridAssembly_sANI95_Comp50Cont5_vs_UHGG_allMAGs"
#magscreen_comp <- "HybridAssembly_sANI98_Comp50Cont5_vs_UHGG_allMAGs"
first_sample <- "ASMDAI005"
final_sample <- "PBSHF003"
#filt <- '*total*strong*relAbund.tsv'
summary_dir <- paste0(counts_dir, prefix, "/summary")
#input <- list.files(summary_dir)[grepl(glob2rx(filt), list.files(summary_dir))]
#input <- "bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund.tsv"
input <- "bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv"
print(input)

# Read in data
counts_table_raw <- read.table(paste0(counts_dir, prefix, "/summary/", input), header = TRUE, sep = "\t")
novel_MAGs <- read.table(paste0(magscreen_dir, magscreen_comp, "/new_species.txt"))$V1
novel_MAGs <- gsub(".fa", "", novel_MAGs)

# Reformat counts table
sample_idx <- match(first_sample, colnames(counts_table_raw))
counts_table <- counts_table_raw[, c(1,2,3, sample_idx:(ncol(counts_table_raw)-2))]
sample_idx2 <- match(first_sample, colnames(counts_table))

# Calculate meta-population prevalence
counts_table$prev_TOTAL <- rowSums(counts_table[, c(sample_idx2:ncol(counts_table))] != 0) / 126

# Calculate population-specific prevalences
header_names <- colnames(counts_table)
pop_names <- c("ASMDAI", "BA0", "BAPDW", "BRUSBG", "MALHUL", "MALRPN", "MALSUL", "PBSHF")
pop_names_grep <- paste(pop_names,"*",sep="")
pop_idx <- list()
for (pop in pop_names_grep) {
  pop_reps <- header_names[grepl(glob2rx(pop), header_names)]
  pop_pref <- gsub("[*]", "", pop)
  pop_idx[[pop_pref]] <- match(pop_reps, header_names)
}
for (pop in pop_names) {
  counts_table[,ncol(counts_table)+1] <- rowSums(counts_table[, c(pop_idx[[pop]])] != 0) / length(pop_idx[[pop]])
  colnames(counts_table)[ncol(counts_table)] <- paste0("prev_", pop)
}

# Calculate standard deviation of population prevalence
prev_list <- colnames(counts_table)[grepl(glob2rx('prev*'), colnames(counts_table))]
prev_list_idx <- match(prev_list, colnames(counts_table))
counts_table$prev_pop_SD <- rowSds(as.matrix(counts_table[,c(prev_list_idx[2:length(prev_list_idx)])]))
counts_table$present_nPops <- rowSums(counts_table[,c(prev_list_idx[2:length(prev_list_idx)])] != 0)

# Partition relative abundance table into novel and non-novel partitions
if(is.na(match("original_bin", colnames(counts_table)))) {
  col_idx2 <- match("Original_bin_ID", colnames(counts_table))
  colnames(counts_table)[col_idx2] <- "original_bin"
}
counts_table_estMAGS <- counts_table[ ! counts_table$original_bin %in% novel_MAGs, ]
#counts_table_estMAGS_downsampled <- counts_table_estMAGS[sample(nrow(counts_table_estMAGS), length(novel_MAGs)), ]
counts_table_novelMAGS <- counts_table[ counts_table$original_bin %in% novel_MAGs, ]
#df_non_novel <- data.frame(partition = "non_novel", mean_rel_abund = relAbund_table_estMAGS$mean_rel_abund)
#df_novel <- data.frame(partition = "novel", mean_rel_abund = relAbund_table_novelMAGS$mean_rel_abund)
#df_partitioned <- rbind(df_non_novel,df_novel)
sample_last_idx <- match(final_sample, colnames(counts_table))
novelMAGS_sumRelAbund <- colSums(counts_table_novelMAGS[,c(sample_idx2:sample_last_idx)])
estMAGS_sumRelAbund <- colSums(counts_table_estMAGS[,c(sample_idx2:sample_last_idx)])
summary_table <- t(cbind(estMAGS_sumRelAbund, novelMAGS_sumRelAbund))
rownames(summary_table) <- c("non-novel MAGs", "novel MAGs")

# Partition prevalence into novel and non-novel partition
df_prev_non_novel <- data.frame(partition = "non_novel", mean_prev = counts_table_estMAGS$prev_TOTAL)
df_prev_novel <- data.frame(partition = "novel", mean_prev = counts_table_novelMAGS$prev_TOTAL)
df_prev_partitioned <- rbind(df_prev_non_novel,df_prev_novel)

# Boxplot
mean_values_prev <- data.frame(group = c("novel", "non-novel"), mean = c(mean(counts_table_novelMAGS$prev_TOTAL), mean(counts_table_estMAGS$prev_TOTAL)))
mean_values_prev$group <- as.factor(mean_values_prev$group)
par(mar = c(5, 5, 3, 3))
#boxplot(df_prev_partitioned$mean_prev ~ df_prev_partitioned$partition, horizontal = TRUE, las=2, col=rgb(0.5,0.65,0.85,0.85), pars=list(outpch = 21, outcol="grey10", outbg = "brown1", outcex=1), cex.axis = 0.5,  xlab="prevalence", ylab="", main="")
#boxplot(df_prev_partitioned$mean_prev ~ df_prev_partitioned$partition, horizontal = TRUE, las=2, col=rgb(0.5,0.65,0.85,0.85), pars=list(outpch = 21, outcol="grey10", outbg = "brown1", outcex=1), cex.axis = 0.5,  ylab="prevalence", xlab="", main="MAG prevalence by novelty")
boxplot(df_prev_partitioned$mean_prev ~ df_prev_partitioned$partition, las=1, col=rgb(0.5,0.65,0.85,0.85), pars=list(outpch = 21, outcol="grey10", outbg = "brown1", outcex=1), cex.axis = 0.5,  ylab="prevalence", xlab="", main="MAG prevalence by novelty")
#points(mean_values_prev$group ~ mean_values_prev$mean, bg = "gold", pch = 23, cex=2)
points(mean_values_prev$mean ~ mean_values_prev$group, bg = "gold", pch = 23, cex=2)

# Two-sample hypothesis testing. Let's use Mann-Whitney U / Wilcoxon rank-sum test, since we don't assume normal distrubution (empirically doesn't follow normal distribution)
wilcox.test(df_prev_non_novel$mean_prev, df_prev_novel$mean_prev, alternative = "two.sided")

# Make summary dataframe
summary_df <- as.data.frame(t(summary_table))
colnames(summary_df) <- c("known_MAGs", "novel_MAGs")
summary_df$pop <- str_split_fixed(rownames(summary_df), "0", 2)[,1]
summary_df$pop <- gsub("3B1B2B3i", "", summary_df$pop)
summary_df$pop <- gsub("3B1B2B3", "", summary_df$pop)
summary_df$pop <- gsub("SINT", "", summary_df$pop)
summary_df$pop <- gsub("BA$", "BA0", summary_df$pop)
summary_df$pop <- as.factor(summary_df$pop)

# Population counts
pop_counts <- table(summary_df$pop)

# Calculate populations means
pop_means <- list()
for (pop in names(pop_counts)) {
  pop_df <- summary_df[summary_df$pop == pop,]
  pop_mean <- mean(pop_df$novel_MAGs)
  pop_means[[pop]] <- pop_mean
}
mean_values <- as.data.frame(unlist(pop_means))
mean_values$pop <- rownames(mean_values)
colnames(mean_values) <- c("mean", "pop")
mean_values$pop <- as.factor(mean_values$pop)

# Plot sum relative abundance of novel vs non-novel MAGS by individual and by population (stacked barplot)
par(mar = c(5, 8, 3, 8))
barplot(summary_table, col=c("sienna1", "skyblue1"), border="white", space=0.1, horiz=TRUE, las=1, cex.names = 0.75, main="Relative abundance of novel vs non-novel MAGs", xlab="Relative abundance", legend.text=TRUE, args.legend=list(x=ceiling(max(summary_table))+0.3, y=ncol(summary_table), bty = "n"))

# Plot boxplots
par(mar = c(5, 5, 3, 3))
boxplot(summary_df$novel_MAGs ~ summary_df$pop, horizontal = TRUE, las=2, col=rgb(0.5,0.65,0.85,0.85), pars=list(outpch = 21, outcol="grey10", outbg = "brown1", outcex=1), cex.axis = 0.5,  xlab="Relative abundance", ylab="", main="Relative abundance of novel MAGs by population")
points(mean_values$pop ~ mean_values$mean, bg = "gold", pch = 23, cex=2)

# Two-sample hypothesis testing. Let's use Mann-Whitney U / Wilcoxon rank-sum test, since we don't assume normal distrubution (empirically doesn't follow normal distribution)
pop_wilcox_df <- pairwise.wilcox.test(summary_df$novel_MAGs, summary_df$pop, p.adjust.method = "fdr", alternative = "two.sided")
pop_wilcox_summary <- pop_wilcox_df$p.value
pop_wilcox_summary # no signficant difference between populations!

# Check if there's correlation between population sample size and relative abundance of novel MAGs
summary_df$pop_n <- summary_df$pop
for (pop in names(pop_counts)) {
  summary_df$pop_n <- gsub(pop, pop_counts[[pop]], summary_df$pop_n)
}
summary_df$pop_n <- as.numeric(summary_df$pop_n)
#plot(summary_df$pop_n, summary_df$novel_MAGs, pch=16, xlab = "population sample size", ylab = "relative abundance of novel MAGs")
#abline(lm(summary_df$novel_MAGs ~ summary_df$pop_n), col="red", lty = 2)
ggplot(summary_df, aes(x=pop_n, y=novel_MAGs)) +
  theme_bw() + rremove("grid") +
  ggtitle("Correlation of the relative abundance of novel MAGs with population sample size") +
  xlab("Population sample size") + ylab("Relative abundance of novel MAGs") +
  geom_smooth(method="lm") +
  geom_point() +
  stat_regline_equation() +
  stat_cor(aes(label=..rr.label..), label.y=0.04)

# Plot MAG metapopulation prevalence (y-axis) vs within-population prevalence standard deviation (x-axis)
par(mar = c(5, 5, 3, 3))
par(mfrow=c(2,1))
#plot(counts_table_novelMAGS$prev_pop_SD, counts_table_novelMAGS$prev_TOTAL, pch=21, col="black", bg=rgb(0,0,0,0.25), ylab = "Total prevalence", xlab = NA, main = "SD of population prevalence vs total (meta-population prevalence)")
#plot(counts_table_estMAGS$prev_pop_SD, counts_table_estMAGS$prev_TOTAL, pch=21, col="black", bg=rgb(0,0,0,0.25), ylab = "Total prevalence", xlab = "SD of population prevalence")
plot(counts_table_novelMAGS$prev_TOTAL, counts_table_novelMAGS$prev_pop_SD, pch=21, col="black", bg=rgb(0,0,0,0.25), ylab = "SD of population prevalence", xlab = NA, main = "Total (meta-population prevalence) vs SD of population prevalence")
plot(counts_table_estMAGS$prev_TOTAL, counts_table_estMAGS$prev_pop_SD, pch=21, col="black", bg=rgb(0,0,0,0.25), ylab = "SD of population prevalence", xlab = "Total prevalence")
ggplot() + theme_bw() +
  geom_point(data = counts_table_estMAGS, aes(x = prev_TOTAL, y = prev_pop_SD), col="grey70") +
  geom_point(data = counts_table_novelMAGS, aes(x = prev_TOTAL, y = prev_pop_SD), col="red") +
  stat_smooth(data = counts_table_novelMAGS, aes(x = prev_TOTAL, y = prev_pop_SD), method = "lm", col = "red",formula = y~poly(x,2),se=F) +
  stat_smooth(data = counts_table_estMAGS, aes(x = prev_TOTAL, y = prev_pop_SD), method = "lm", col = "black",formula = y~poly(x,2),se=F) +
  geom_line()

# Compare population-specificity of novel vs non-novel MAGs
counts_table_novelMAGS_pop <- counts_table_novelMAGS[counts_table_novelMAGS$present_nPops == 1, ]
counts_table_novelMAGS_pop <- counts_table_novelMAGS_pop[,c(131:138)]
counts_table_novelMAGS_pop$pop <- apply(counts_table_novelMAGS_pop, 1, function(x) paste(colnames(counts_table_novelMAGS_pop)[which(x > 0)], collapse = ", "))
counts_table_novelMAGS_pop$pop <- gsub("prev_", "", counts_table_novelMAGS_pop$pop)
counts_table_novelMAGS_popTab <- as.data.frame(table(counts_table_novelMAGS_pop$pop))
counts_table_novelMAGS_popTab[,c("2","3","4","5","6","7","8")] <- 0
rownames(counts_table_novelMAGS_popTab) <- counts_table_novelMAGS_popTab$Var1
counts_table_novelMAGS_popTab <- counts_table_novelMAGS_popTab[,-1]
counts_table_novelMAGS_popTab <- rbind(counts_table_novelMAGS_popTab, table(counts_table_novelMAGS$present_nPops))
rownames(counts_table_novelMAGS_popTab)[nrow(counts_table_novelMAGS_popTab)] <- "total"
colnames(counts_table_novelMAGS_popTab)[1] <- "1"
counts_table_novelMAGS_popTab[nrow(counts_table_novelMAGS_popTab),1] <- 0

#cols <- c(brewer.pal(7, "Pastel2"), "skyblue1")
#cols <- c("#FF6347", "#277B45", "#4CFF8E", "#FDAC10", "#7C378D", "#DF96B4", "#B494D5", "skyblue1")
pop_cols <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Pilot ms/pop_colours.tsv", header = TRUE, sep = "\t", comment.char="") 
cols <- c(pop_cols$col[1:7], "skyblue1")
barplot(table(counts_table_estMAGS$present_nPops), border="white", ylab = "Frequency", main = "Population-specificity of novel (blue) vs non-novel (orange) MAGs", col = "sienna1", space=rep(0,7))
#barplot(table(counts_table_novelMAGS$present_nPops), border="white", ylab = "Frequency", xlab = "Number of populations MAG is present in", col = "skyblue1", space=rep(0,7))
barplot(as.matrix(counts_table_novelMAGS_popTab), border="white", col=cols, ylab = "Frequency", xlab = "Number of populations MAG is present in", space=rep(0,7))

# Perform K-S test. Alternative null hypothesis can be one of "less", "greater", "two.sided"
ks_res <- ks.test(counts_table_estMAGS$present_nPops, counts_table_novelMAGS$present_nPops, alternative = "two.sided", exact = TRUE)

# Write out table of MAG population-specificity
all(colnames(counts_table_estMAGS) == colnames(counts_table_novelMAGS))
counts_table_allMAGs <- rbind(counts_table_estMAGS, counts_table_novelMAGS)
#write.table(counts_table_allMAGs, "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/Novel vs non-novel comparisons/MAG_pop_specificity_revised.txt", sep = "\t", quote = FALSE, row.names = FALSE)
