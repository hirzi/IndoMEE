### This script compares and plots the village-level prevalence of sporulating vs non-sporulating MAGs

# Load libraries
library(gtools)
library(matrixStats)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(stringr)

# Set working directories
counts_dir <- "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assemblies/On_hostDecontaminated_fastqs/"

# Define parameters
#assembly_list <- c("CoAssembly_draft1", "CoAssembly_sANI98_draft2", "HybridAssembly_Comp50Cont5_draft2", "HybridAssembly_sANI98_Comp50Cont5_draft2", "Hybrid_UHGG_Assembly_Comp50Cont5_draft2", "Hybrid_UHGG_SPMP_Assembly_Comp50Cont5_draft2", "IndAssembly_draft2", "IndAssembly_draft2_sANI98", "SPMP", "UHGG")
prefix <- "HybridAssembly_Comp50Cont5_draft2"
first_sample <- "ASMDAI005"
final_sample <- "PBSHF003"
#filt <- '*total*strong*relAbund.tsv'
summary_dir <- paste0(counts_dir, prefix, "/summary")
#input <- list.files(summary_dir)[grepl(glob2rx(filt), list.files(summary_dir))]
input <- "bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv"
print(input)

# Read in data
counts_table_raw <- read.table(paste0(counts_dir, prefix, "/summary/", input), header = TRUE, sep = "\t")

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

# Define run parameters
filter_MAGs <- TRUE
MAG_type_thres <- "binary" # counts or binary

# Filter by MAGs (e.g. sporing / non-sporing)
counts_table_filt <- counts_table
if(filter_MAGs==TRUE) {
  if(MAG_type_thres=="binary"){
    # Get spo0A MAGs. First, run functionalEnrichmentAnalysis.R
    feat.spo0A <- feat.all[grep("spo0A", rownames(feat.all)), ]
    colnames(feat.spo0A) <- gsub("treatment_", "", colnames(feat.spo0A))
    colnames(feat.spo0A) <- gsub("control_", "", colnames(feat.spo0A))
    feat.spo0A <- as.data.frame(t(feat.spo0A))
    colnames(feat.spo0A) <- "spo0A"
    feat.spo0A$original_bin <- rownames(feat.spo0A)
    #spore_MAGs <- feat.spo0A[feat.spo0A$spo0A==1,]$original_bin
    spore_MAGs <- feat.spo0A[feat.spo0A$spo0A>=1,]$original_bin
    nonSpore_MAGs <- feat.spo0A[feat.spo0A$spo0A==0,]$original_bin
    #spore_MAGs <- sample(spore_MAGs, min(length(spore_MAGs), length(nonSpore_MAGs)))
    #nonSpore_MAGs <- sample(nonSpore_MAGs, min(length(spore_MAGs), length(nonSpore_MAGs)))
  } else if(MAG_type_thres=="counts"){
    # Alternatively, partition by spore counts. Get spore counts for Bacilotta A. First, run functionalEnrichmentAnalysis.R, followed by Sporulation_features_presenceAbsence_logisticRegression.R (till the appropriate row).
    # Consider partition e.g. top 200 vs bottom 200.
    p <- phyla[3]
    print(p)
    phylum_df <- phyla_df[phyla_df$phylum == p, ]
    feat.sel.prev.phylum <- merge(phylum_df, feat.sel.prev, by = "original_bin")
    feat.sel.sum <- feat.sel.prev.phylum[,c(2:(ncol(feat.sel.prev.phylum)-1))]
    feat.sel.sum <- feat.sel.sum[,-c(1)]
    feat.sel.sum$sum <- rowSums(feat.sel.sum)
    feat.sel.sum$prevalence <-feat.sel.prev.phylum$prevalence
    feat.sum.BacA <- as.data.frame(feat.sel.sum$sum)
    rownames(feat.sum.BacA) <- feat.sel.prev.phylum$original_bin
    colnames(feat.sum.BacA) <- c("sum")
    feat.sum.BacA$original_bin <- rownames(feat.sum.BacA)
    spore_MAGs <- feat.sum.BacA[feat.sum.BacA$sum>23,]$original_bin
    nonSpore_MAGs <- feat.sum.BacA[feat.sum.BacA$sum<10,]$original_bin
    #spore_MAGs <- sample(spore_MAGs, min(length(spore_MAGs), length(nonSpore_MAGs)))
    #nonSpore_MAGs <- sample(nonSpore_MAGs, min(length(spore_MAGs), length(nonSpore_MAGs)))
  }
  # And filter
  rownames(counts_table_filt) <- counts_table_filt$original_bin
  counts_table_sporers <- counts_table_filt[rownames(counts_table_filt) %in% spore_MAGs, ]
  counts_table_nonsporers <- counts_table_filt[rownames(counts_table_filt) %in% nonSpore_MAGs, ]
}

# Compare population-specificity of sporulating vs non-sporulating MAGs
counts_table_nonsporers_pop <- counts_table_nonsporers[counts_table_nonsporers$present_nPops == 1, ]
counts_table_nonsporers_pop <- counts_table_nonsporers_pop[,c(131:138)]
counts_table_nonsporers_pop$pop <- apply(counts_table_nonsporers_pop, 1, function(x) paste(colnames(counts_table_nonsporers_pop)[which(x > 0)], collapse = ", "))
counts_table_nonsporers_pop$pop <- gsub("prev_", "", counts_table_nonsporers_pop$pop)
counts_table_nonsporers_popTab <- as.data.frame(table(counts_table_nonsporers_pop$pop))
counts_table_nonsporers_popTab[,c("2","3","4","5","6","7","8")] <- 0
rownames(counts_table_nonsporers_popTab) <- counts_table_nonsporers_popTab$Var1
counts_table_nonsporers_popTab <- counts_table_nonsporers_popTab[,-1]
counts_table_nonsporers_popTab <- rbind(counts_table_nonsporers_popTab, table(counts_table_nonsporers$present_nPops))
rownames(counts_table_nonsporers_popTab)[nrow(counts_table_nonsporers_popTab)] <- "total"
colnames(counts_table_nonsporers_popTab)[1] <- "1"
counts_table_nonsporers_popTab[nrow(counts_table_nonsporers_popTab),1] <- 0

#cols <- c(brewer.pal(7, "Pastel2"), "skyblue1")
#cols <- c("#FF6347", "#277B45", "#4CFF8E", "#FDAC10", "#7C378D", "#DF96B4", "#B494D5", "skyblue1")
pop_cols <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Pilot ms/pop_colours.tsv", header = TRUE, sep = "\t", comment.char="") 
#cols <- c(pop_cols$col[1:7], "skyblue1")
cols <- c(pop_cols$col, "skyblue1")
par(mar = c(5, 5, 3, 3))
par(mfrow=c(2,1))
barplot(table(counts_table_sporers$present_nPops), border="white", ylab = "Frequency", ylim = c(0,120), main = "Population-specificity of sporing (orange) vs non-sporing (blue) MAGs", col = "sienna1", space=rep(0,7))
#barplot(table(counts_table_nonsporers$present_nPops), border="white", ylab = "Frequency", xlab = "Number of populations MAG is present in", col = "skyblue1", space=rep(0,7))
barplot(as.matrix(counts_table_nonsporers_popTab), border="white", col=cols, ylab = "Frequency", ylim = c(0,120), xlab = "Number of populations MAG is present in", space=rep(0,7))
# Side-by-side barplot
counts_table_sideByside <- rbind(table(counts_table_sporers$present_nPops), table(counts_table_nonsporers$present_nPops))
rownames(counts_table_sideByside) <- c("sporing", "non-sporing")
par(mfrow=c(1,1))
barplot(counts_table_sideByside, beside=T, legend.text=T, col=c("#afafaf" , "skyblue1"), ylab="Frequency", space=c(0,0.7))

# Perform K-S test. Alternative null hypothesis can be one of "less", "greater", "two.sided"
ks_res <- ks.test(counts_table_sporers$present_nPops, counts_table_nonsporers$present_nPops, alternative = "two.sided", exact = TRUE)
