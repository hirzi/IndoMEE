# This scripts compares the lengths and gene counts of novel vs known MAGs
# It can compare raw, filtered (by completeness) and normalised (by completeness) metrics.
# Current results first applies a 90% completeness hard filter; however we note that the distribution of completeness is still significantly skewed in favour of known MAGs
# Thus, we additionally normalise by completeness. We don't normalise all MAGs (i.e. without prior hard filter) because completion metrics are approximate, and we use it primarily to retain only the highest quality, most confident/complete genomes.

# Define paths
feat.path = "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/Functional enrichment/HybridAssembly_Comp50Cont5_Summary/"
descr.path <- "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/Functional enrichment/"
magscreen.path <- "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/magscreen/Novel_MAGs/"
all.MAGs <- scan("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/Functional enrichment/dereplicated_genomes.txt", what = "")
MAG.metadata <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Metadata/metawrap_HybridAssembly_Comp50Cont5_draft2_metadata.tsv", header=TRUE)
file <- paste0(feat.path, "eggnog_emapper_summary.tsv")
feat.data = read.delim(file, header=FALSE, stringsAsFactors = FALSE)[,c(1,10)][-1, ]
feat.data = read.delim(file, header=TRUE, stringsAsFactors = FALSE)[,c(1,7)]

# Load data - MAGs
magscreen_comp <- "HybridAssembly_sANI95_Comp50Cont5_vs_UHGG_allMAGs"
novel.MAGs <- scan(paste0(magscreen.path, magscreen_comp, "/new_species.txt"), what = "")
novel.MAGs <- gsub(".fa", "", novel.MAGs)
all.MAGs <- gsub(".fa", "", all.MAGs)
known.MAGs <- all.MAGs[which(!all.MAGs %in% novel.MAGs)]
#known.MAGs <- sample(known.MAGs, length(novel.MAGs), replace = FALSE)
control.MAGs <- known.MAGs
treatment.MAGs <- novel.MAGs

# Filter to MAGs above certain completeness
compl_filter <- TRUE
compl_theshold <- 90
if (compl_filter == TRUE) {
  MAG.metadata.filt <- MAG.metadata[MAG.metadata$completeness > compl_theshold, ]
  all.MAGs <- all.MAGs[which(all.MAGs %in% MAG.metadata.filt$original_bin)]
  treatment.MAGs <- treatment.MAGs[which(treatment.MAGs %in% MAG.metadata.filt$original_bin)]
  control.MAGs <- control.MAGs[which(control.MAGs %in% MAG.metadata.filt$original_bin)]
}

# Extract features
proteins.prefixes <-  sub("_[^_]*$", "", feat.data[,1])
# For treatment MAGs
i <- 1
for (MAG in treatment.MAGs) {
  feat.treatment.MAG = as.vector(as.matrix(feat.data[which(proteins.prefixes %in% MAG),-1]))
  feat.treatment.MAG = feat.treatment.MAG[which(feat.treatment.MAG != "")]
  multi_feats <- feat.treatment.MAG[nchar(feat.treatment.MAG) > 1]
  single_feats <- feat.treatment.MAG[nchar(feat.treatment.MAG) == 1]
  multi_feat.unlist <- unlist(sapply(multi_feats, function(x) unlist(strsplit(x, split = ""))))
  feat.treatment.MAG <- c(single_feats, multi_feat.unlist)
  if (length(feat.treatment.MAG) > 0) {
    feat.treatment.MAG.df <- data.frame(table(feat.treatment.MAG))
  } else if (length(feat.treatment.MAG) == 0) {
    feat.treatment.MAG.df <- data.frame(feature=factor(),Freq=integer())
  }
  colnames(feat.treatment.MAG.df) = c("feature", paste0("treatment_",MAG))
  if (i == 1) {
    feat.treatment.df <- feat.treatment.MAG.df
  } else {
    feat.treatment.df = merge(feat.treatment.df, feat.treatment.MAG.df, by="feature", all = TRUE)
  }
  i <- i + 1
}
feat.treatment.df[is.na(feat.treatment.df)] <- 0

# For control MAGs
i <- 1
for (MAG in control.MAGs) {
  feat.control.MAG = as.vector(as.matrix(feat.data[which(proteins.prefixes %in% MAG),-1]))
  feat.control.MAG = feat.control.MAG[which(feat.control.MAG != "")]
  multi_feats <- feat.control.MAG[nchar(feat.control.MAG) > 1]
  single_feats <- feat.control.MAG[nchar(feat.control.MAG) == 1]
  multi_feat.unlist <- unlist(sapply(multi_feats, function(x) unlist(strsplit(x, split = ""))))
  feat.control.MAG <- c(single_feats, multi_feat.unlist)
  if (length(feat.control.MAG) > 0) {
    feat.control.MAG.df = data.frame(table(feat.control.MAG))
  } else if (length(feat.control.MAG) == 0) {
    feat.control.MAG.df <- data.frame(feature=factor(),Freq=integer())
  }
  colnames(feat.control.MAG.df) = c("feature", paste0("control_",MAG))
  if (i == 1) {
    feat.control.df <- feat.control.MAG.df
  } else {
    feat.control.df = merge(feat.control.df, feat.control.MAG.df, by="feature", all = TRUE)
  }
  i <- i + 1
}
feat.control.df[is.na(feat.control.df)] <- 0

# Merge and format data
feat.all = merge(feat.treatment.df, feat.control.df, by="feature", all = TRUE)
feat.all[is.na(feat.all)] <- 0
rownames(feat.all) <- feat.all$feature
feat.all <- feat.all[,-1]
treatment_idx <- which(colnames(feat.all) %in% paste0("treatment_", treatment.MAGs))
control_idx <- which(colnames(feat.all) %in% paste0("control_", control.MAGs))
#rm(feat.data, feat.control.df, feat.control.MAG.df, feat.control.MAG, feat.treatment.df, feat.treatment.MAG.df, feat.treatment.MAG)

# Get total gene counts and metadata per MAG
geneCounts <- data.frame(colSums(feat.all))
colnames(geneCounts)[1] <- "gene_counts"
geneCounts$original_bin <- str_split_fixed(rownames(geneCounts), "_", 2)[,2]
geneCounts$idx  <- 1:nrow(geneCounts)
geneCounts <- merge(geneCounts, MAG.metadata, by = "original_bin")
geneCounts <- geneCounts[,c(1,2,3,24,25,26)]
geneCounts <- geneCounts[order(geneCounts$idx), ]
# Check ordering is correct
all(geneCounts$original_bin == (str_split_fixed(colnames(feat.all), "_", 2)[,2]))

# Check whether distribution of completeness, post length filter, is significantly different between novel and known MAGs, which if so, suggest need for completeness normalisation.
treatment_geneCounts <- geneCounts[treatment_idx, ]
control_geneCounts <- geneCounts[control_idx, ]
compl_comp90 <- rbind(data.frame(compl=control_geneCounts$completeness, group='known'),
                      data.frame(compl=treatment_geneCounts$completeness, group='novel'))
boxplot(compl_comp90$compl ~ compl_comp90$group, las=1, col=rgb(0.5,0.65,0.85,0.85), pars=list(outpch = 21, outcol="grey10", outbg = "brown1", outcex=1), cex.axis = 0.5,  ylab="length (bp)", xlab="", main="MAG completeness")
wilcox.test(control_geneCounts$completeness, treatment_geneCounts$completeness, alternative = "greater")

# Normalise by completion rate
geneCounts$gene_counts_norm <- geneCounts$gene_counts / (geneCounts$completeness / 100)
geneCounts$length_norm <- geneCounts$length / (geneCounts$completeness / 100)

# Partition into known and novel
treatment_geneCounts <- geneCounts[treatment_idx, ]
control_geneCounts <- geneCounts[control_idx, ]

### GENE COUNTS ###
# Output means
mean(treatment_geneCounts$gene_counts)
mean(control_geneCounts$gene_counts)

mean(treatment_geneCounts$gene_counts_norm)
mean(control_geneCounts$gene_counts_norm)

# Plot boxplots
gene_counts_bp <- rbind(data.frame(counts=control_geneCounts$gene_counts, group='known'),
                    data.frame(counts=treatment_geneCounts$gene_counts, group='novel'))
gene_counts_norm_bp <- rbind(data.frame(counts=control_geneCounts$gene_counts_norm, group='known'),
                     data.frame(counts=treatment_geneCounts$gene_counts_norm, group='novel'))

boxplot(gene_counts_bp$counts ~ gene_counts_bp$group, las=1, col=rgb(0.5,0.65,0.85,0.85), pars=list(outpch = 21, outcol="grey10", outbg = "brown1", outcex=1), cex.axis = 0.5,  ylab="gene count", xlab="", main="MAG gene counts")
boxplot(gene_counts_norm_bp$counts ~ gene_counts_norm_bp$group, las=1, col=rgb(0.5,0.65,0.85,0.85), pars=list(outpch = 21, outcol="grey10", outbg = "brown1", outcex=1), cex.axis = 0.5,  ylab="gene count", xlab="", main="MAG gene counts")

# Two-sample hypothesis testing. Let's use Mann-Whitney U / Wilcoxon rank-sum test, since we don't assume normal distribution
wilcox.test(control_geneCounts$gene_counts, treatment_geneCounts$gene_counts, alternative = "greater")
wilcox.test(control_geneCounts$gene_counts_norm, treatment_geneCounts$gene_counts_norm, alternative = "greater")

### LENGTH ###
# Output means
mean(treatment_geneCounts$length)
mean(control_geneCounts$length)

mean(treatment_geneCounts$length_norm)
mean(control_geneCounts$length_norm)

# Plot boxplots
length <- rbind(data.frame(length=control_geneCounts$length, group='known'),
                       data.frame(length=treatment_geneCounts$length, group='novel'))
length_norm <- rbind(data.frame(length=control_geneCounts$length_norm, group='known'),
                            data.frame(length=treatment_geneCounts$length_norm, group='novel'))

boxplot(length$length ~ length$group, las=1, col=rgb(0.5,0.65,0.85,0.85), pars=list(outpch = 21, outcol="grey10", outbg = "brown1", outcex=1), cex.axis = 0.5,  ylab="length (bp)", xlab="", main="MAG length")
boxplot(length_norm$length ~ length_norm$group, las=1, col=rgb(0.5,0.65,0.85,0.85), pars=list(outpch = 21, outcol="grey10", outbg = "brown1", outcex=1), cex.axis = 0.5,  ylab="length (bp)", xlab="", main="MAG length")

# Two-sample hypothesis testing. Let's use Mann-Whitney U / Wilcoxon rank-sum test, since we don't assume normal distribution
wilcox.test(control_geneCounts$length, treatment_geneCounts$length, alternative = "greater")
wilcox.test(control_geneCounts$length_norm, treatment_geneCounts$length_norm, alternative = "greater")

