# This scripts plots regressions of MAG (significantly enriched) sporulation features against MAG prevalence
# Run functionalEnrichmentAnalysis.R before running this script.
# It first performs logistic regressions on the presence-absence of (significantly-associated) sporulation features against prevalence, per feature.
# It then performs a standard regression on the count of (presences) of (significantly-associated) sporulation features against prevalence.

# Load libraries
library(ggplot2)
library(stringr)
library(gtools)

# Indicate whether to plot
plot_all <- TRUE

# Read in MAG metadata
MAG.metadata2 <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Phylogenetics/HybridAssembly_Comp50Cont5_draft2/Metadata_HybridAssembly_Comp50Cont5_draft2_revised25Nov_MaaslinNeedsToBeUpdated.txt", comment.char="", sep = "\t", header = TRUE)

# Get MAG presence-absence of significantly associated sporulation-associated features
feat.all.temp <- feat.all
feat.all.temp$feature <- rownames(feat.all.temp)
feat.sel <- merge(res.sig, feat.all.temp, by="feature")
if (feat == "pfam") {
  feat.sel <- feat.sel[grep("Spo", feat.sel$feature), ]
  rownames(feat.sel) <- feat.sel$feature
  feat.sel <- feat.sel[,-c(1:5)]
} else if (feat == "GO") {
  feat.sel1 <- feat.sel[grep("spore", feat.sel$name), ]
  feat.sel2 <- feat.sel[grep("sporu", feat.sel$name), ]
  feat.sel <- merge(feat.sel1, feat.sel2, all=TRUE)
  rm(feat.sel1, feat.sel2)
  feat.sel$feature <- gsub("00", "GO_00", feat.sel$feature)
  rownames(feat.sel) <- feat.sel$feature
  feat.sel <- feat.sel[,-c(1:8)]
} else if (feat == "kegg") {
  feat.sel$feature <- str_split_fixed(feat.sel$feature, ";", n = Inf)[,1]
  feat.sel1 <- feat.sel[grep("spore", feat.sel$description), ]
  feat.sel2 <- feat.sel[grep("sporu", feat.sel$description), ]
  feat.sel <- merge(feat.sel1, feat.sel2, all=TRUE)
  rm(feat.sel1, feat.sel2)
  rownames(feat.sel) <- feat.sel$feature
  feat.sel <- feat.sel[,-c(1:5)]
}
feat.sel.t <- as.data.frame(t(feat.sel))
rownames(feat.sel.t) <- gsub("treatment_", "", rownames(feat.sel.t))
rownames(feat.sel.t) <- gsub("control_", "", rownames(feat.sel.t))
feat.sel.t$original_bin <- rownames(feat.sel.t)
prev.df <- MAG.metadata2[,c("original_bin", "prevalence")]
feat.sel.prev <- merge(feat.sel.t, prev.df, by = "original_bin")
feat.sel.prev[,c(2:(ncol(feat.sel.prev)-1))] <- sapply(feat.sel.prev[,c(2:(ncol(feat.sel.prev)-1))], as.integer)

# Output sporulation gene counts table
feat.sel.alt <- feat.sel.t[,-c(ncol(feat.sel.t))]
feat.sel.alt$n_spore_genes <- rowSums(feat.sel.alt)
#write.table(feat.sel.alt, paste0("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/Functional enrichment/", feat, "_", feat_type, "_sporeGeneCounts.txt"), quote = FALSE, row.names = TRUE, sep = "\t")

# Plot logistic regressions of MAG sporulation feature presence-absences vs prevalence.
if (plot_all == TRUE) {
  pdf(file = paste0("/Users/hl636/Desktop/", "sporePresenceAbsence_", feat, "_", feat_type, "_byFeature.pdf"), width = 8, height = 5)
  for (i in seq(2,(ncol(feat.sel.prev)-1))){
    print(i)
    resp_var <- colnames(feat.sel.prev)[i]
    log_form <- as.formula(paste(resp_var, "prevalence", sep = " ~ "))
    # Fit logistic regression model
    model <- glm(log_form, data=feat.sel.prev, family=binomial)
    # Define new data frame that contains predictor variable
    newdata <- data.frame(prevalence=seq(min(feat.sel.prev$prevalence), max(feat.sel.prev$prevalence),len=500))
    # Use fitted model to predict values of response variable
    newdata$resp_var = predict(model, newdata=newdata, type="response")
    # Plot logistic regression curve
    plot(log_form, data=feat.sel.prev, col="steelblue", ylim=c(0,1), xlab="prevalence", ylab="0:absence - 1:presence", main=paste0("All MAGs: ", colnames(feat.sel.prev)[i]))
    lines(resp_var ~ prevalence, newdata, lwd=2)
  }
  dev.off()
}


############## Across all phyla ##############
# Read in number of populations MAG is present in
counts_table_allMAGs <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/Novel vs non-novel comparisons/MAG_pop_specificity_revised.txt", sep = "\t", header = TRUE)

# Plot the regression of the count (i.e. number of presences) of (significantly-associated) sporulation features against prevalence.
feat.sel.prev <- merge(counts_table_allMAGs, feat.sel.prev, by = "original_bin")
idx1 <- match("present_nPops", colnames(feat.sel.prev))
feat.sel.prev <- feat.sel.prev[, c(1, (idx1+1):ncol(feat.sel.prev), idx1)]
feat.sel.sum <- feat.sel.prev[,c(2:(ncol(feat.sel.prev)-2))]
feat.sel.sum$sum <- rowSums(feat.sel.sum)
feat.sel.sum$prevalence <-feat.sel.prev$prevalence
feat.sel.sum$present_nPops <-feat.sel.prev$present_nPops
# And plot
if (plot_all == TRUE) {
  pdf(file = paste0("/Users/hl636/Desktop/", "sporeCounts_", feat, "_", feat_type, ".pdf"), width = 3.3, height = 3)
  print(ggscatter(feat.sel.sum, x = "sum", y = "prevalence", add = "reg.line", size = 1, color = "black", fill="grey60", shape = 21, add.params = list(color = "red2", linetype= "dashed")) + theme_bw() +
          stat_cor(label.x = 0.6, label.y = 0.1) + ylim(0,1) + xlim(min(feat.sel.sum$sum), max(feat.sel.sum$sum)) + ggtitle("All MAGs: prevalence - sporulation gene counts") + xlab("counts") +
          stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), label.x = 0.6, label.y = 0.05) +
          rremove("grid"))
  print(ggscatter(feat.sel.sum, x = "sum", y = "present_nPops", add = "reg.line", size = 1, color = rgb(0.2,0.2,0.2,0.2), fill=rgb(0.2,0.2,0.2,0.2), shape = 21, add.params = list(color = "red2", linetype= "dashed")) + theme_bw() +
          stat_cor(label.x = min(feat.sel.sum$sum), label.y = min(feat.sel.sum$present_nPops) + 2) + ylim(1,8) + ggtitle("All MAGs: present_nPops - sporulation gene counts") + ylab("present_nPops") + xlab("counts") +
          stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), label.x = min(feat.sel.sum$sum), label.y = min(feat.sel.sum$present_nPops) + 1) +
          rremove("grid"))
  dev.off()
}


############## Run by phylum ##############

prefix <- "HybridAssembly_Comp50Cont5_draft2"
metadata_df <- read.table(paste0("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Phylogenetics/HybridAssembly_Comp50Cont5_draft2/Metadata_", prefix, "_revised25Nov_MaaslinNeedsToBeUpdated.txt"), quote = "", comment.char = "", sep = "\t", header = TRUE)
phyla_df <- metadata_df[,c("original_bin", "phylum")]
phyla <- mixedsort(unique(phyla_df$phylum))

# Plot logistic regressions of MAG sporulation feature presence-absences vs prevalence.
if (plot_all == TRUE) {
  pdf(file = paste0("/Users/hl636/Desktop/", "sporePresenceAbsence_", feat, "_", feat_type, "_byPhyla_byFeature.pdf"), width = 8, height = 5)
  par(mar = c(4, 4, 3, 3), mfrow = c(1, 1))
  for (p in phyla) {
    print(p)
    phylum_df <- phyla_df[phyla_df$phylum == p, ]
    feat.sel.prev.phylum <- merge(phylum_df, feat.sel.prev, by = "original_bin")
    for (i in seq(3,(ncol(feat.sel.prev.phylum)-2))){
      print(i)
      resp_var <- colnames(feat.sel.prev.phylum)[i]
      log_form <- as.formula(paste(resp_var, "prevalence", sep = " ~ "))
      # Fit logistic regression model
      model <- glm(log_form, data=feat.sel.prev.phylum, family=binomial)
      # Define new data frame that contains predictor variable
      newdata <- data.frame(prevalence=seq(min(feat.sel.prev.phylum$prevalence), max(feat.sel.prev.phylum$prevalence),len=500))
      # Use fitted model to predict values of response variable
      newdata$resp_var = predict(model, newdata=newdata, type="response")
      # Plot logistic regression curve
      plot(log_form, data=feat.sel.prev.phylum, col="steelblue", main=paste0(p, ": ", colnames(feat.sel.prev.phylum)[i]), ylim=c(0,1), xlab="prevalence")
      lines(resp_var ~ prevalence, newdata, lwd=2)
    }
  }
  dev.off()
}

# Plot the regression of the count (i.e. number of presences) of (significantly-associated) sporulation features against prevalence.
# Print significant results. Scale 180 in AI%
# Interesting ones are Bacilotta*, Bacilotta_A*, Bacilotta_B, Bacilotta_C, Bacteroidota*, Spirochaetota*; only those denoted by * are significant
#phyla_sel <- phyla[c(2,3,4,5,6,16)]
if (plot_all == TRUE) {
  pdf(file = paste0("/Users/hl636/Desktop/", "sporeCounts_", feat, "_", feat_type, "_byPhyla_final.pdf"), width = 3.3, height = 3)
  for (p in phyla) {
    print(p)
    phylum_df <- phyla_df[phyla_df$phylum == p, ]
    feat.sel.prev.phylum <- merge(phylum_df, feat.sel.prev, by = "original_bin")
    feat.sel.sum <- feat.sel.prev.phylum[,c(3:(ncol(feat.sel.prev.phylum)-2))]
    feat.sel.sum$sum <- rowSums(feat.sel.sum)
    feat.sel.sum$prevalence <-feat.sel.prev.phylum$prevalence
    feat.sel.sum$present_nPops <-feat.sel.prev.phylum$present_nPops
    # And plot
    print(ggscatter(feat.sel.sum, x = "sum", y = "prevalence", add = "reg.line", size = 1, color = "black", fill="grey60", shape = 21, add.params = list(color = "red2", linetype= "dashed")) + theme_bw() +
            stat_cor(label.x = min(feat.sel.sum$sum), label.y = 0.1) + ylim(min(feat.sel.sum$prevalence), max(feat.sel.sum$prevalence)) + xlim(min(feat.sel.sum$sum), max(feat.sel.sum$sum)) + ggtitle(paste0(p, ": sporulation gene counts")) + xlab("counts") +
            stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), label.x = min(feat.sel.sum$sum), label.y = 0.05) +
            rremove("grid"))
    print(ggscatter(feat.sel.sum, x = "sum", y = "present_nPops", add = "reg.line", size = 1, color = rgb(0.2,0.2,0.2,0.1), fill=rgb(0.2,0.2,0.2,0.1), shape = 21, add.params = list(color = "red2", linetype= "dashed")) + theme_bw() +
            stat_cor(label.x = min(feat.sel.sum$sum), label.y = 2) + ylim(1, 8) + xlim(min(feat.sel.sum$sum), max(max(feat.sel.sum$sum),1)) + ggtitle(paste0(p, ": sporulation gene counts")) + ylab("present_nPops") +
            stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), label.x = min(feat.sel.sum$sum), label.y = 1) +
            rremove("grid"))
  }
  dev.off()
}

### Multiple comparison adjustment on Bacillota phyla p-values
# Prevalence
p_vals_bacillota <- c(0.00072, 0.13, 0.45, 3.4e-5) # KEGG; Bacilotta_A, Bacilotta_B, Bacilotta_C, Bacilotta
p_vals_allPhyla <- c(0.19, 3.4e-5, 0.00072, 0.13, 0.45, 0.013, NA, 0.52, NA, NA, NA, 0.5, NA, NA, 0.056, 0.022, NA, NA, NA) # ordered alphabetically
p_vals_bacillota <- c(1.8e-08, 0.5, 0.75, 0.00015) # Pfam; Bacilotta_A, Bacilotta_B, Bacilotta_C, Bacilotta
p.adjust(p=p_vals_bacillota, method = "fdr", n=length(p_vals_bacillota))
p.adjust(p=p_vals_allPhyla, method = "fdr", n=length(p_vals_allPhyla))
# Number of populations present in
p_vals_bacillota <- c(0.00019, 0.2, 0.14, 0.046) # Bacilotta_A, Bacilotta_B, Bacilotta_C, Bacilotta ss
p_vals_allPhyla <- c(0.13, 0.046, 0.00019, 0.2, 0.14, 0.0026, NA, 0.41, NA, NA, NA, 0.72, NA, NA, 0.009, 0.01, NA, NA, NA) # ordered alphabetically
p.adjust(p=p_vals_bacillota, method = "fdr", n=length(p_vals_bacillota))
p.adjust(p=p_vals_allPhyla, method = "fdr", n=length(p_vals_allPhyla))

############## Run pan Bacillota ##############
# Plot the regression of the count (i.e. number of presences) of (significantly-associated) sporulation features against prevalence.
# For pan-Bacilotta (i.e. Bacillota + Bacillota_A + Bacillota_B + Bacillota_C)
phylum_df <- phyla_df[grep("Bacillota", phyla_df$phylum), ]
feat.sel.prev.phylum <- merge(phylum_df, feat.sel.prev, by = "original_bin")
feat.sel.sum <- feat.sel.prev.phylum[,c(3:(ncol(feat.sel.prev.phylum)-2))]
feat.sel.sum$sum <- rowSums(feat.sel.sum)
feat.sel.sum$prevalence <-feat.sel.prev.phylum$prevalence
feat.sel.sum$present_nPops <-feat.sel.prev.phylum$present_nPops
# And plot
if (plot_all == TRUE) {
  pdf(file = paste0("/Users/hl636/Desktop/", "sporeCounts_", feat, "_", feat_type, "_panBacillota_significant_final.pdf"), width = 3.3, height = 3)
  print(ggscatter(feat.sel.sum, x = "sum", y = "prevalence", add = "reg.line", size = 1, color = "black", fill="grey60", shape = 21, add.params = list(color = "red2", linetype= "dashed")) + theme_bw() +
          stat_cor(label.x = 1, label.y = 0.8) + xlim(min(feat.sel.sum$sum), max(feat.sel.sum$sum)) + ylim(min(feat.sel.sum$prevalence), max(feat.sel.sum$prevalence)) + ggtitle("Bacillota: sporulation gene counts") + ylab("prevalence") +
          stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), label.x = 1, label.y = 0.6) +
          rremove("grid"))
  print(ggscatter(feat.sel.sum, x = "sum", y = "present_nPops", add = "reg.line", size = 1, color = rgb(0.2,0.2,0.2,0.2), fill=rgb(0.2,0.2,0.2,0.2), shape = 21, add.params = list(color = "red2", linetype= "dashed")) + theme_bw() +
          stat_cor(label.x = min(feat.sel.sum$sum), label.y = min(feat.sel.sum$present_nPops) + 2) + ylim(1,8) + ggtitle("Bacillota: sporulation gene counts") + ylab("present_nPops") +
          stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), label.x = min(feat.sel.sum$sum), label.y = min(feat.sel.sum$present_nPops) + 1) +
          rremove("grid"))
  dev.off()
}

############## Run other clades ##############
# Plot the regression of the count (i.e. number of presences) of (significantly-associated) sporulation features against prevalence.
# For Bacilotta spore-forming clade
bacillota_spore_clade <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/Functional enrichment/Bacillota_sporeFormingClade.txt", header = FALSE)
colnames(bacillota_spore_clade) <- "original_bin"
feat.sel.prev.bacillotaSporeClade <- merge(bacillota_spore_clade, feat.sel.prev, by = "original_bin")
feat.sel.sum <- feat.sel.prev.bacillotaSporeClade[,c(2:(ncol(feat.sel.prev.bacillotaSporeClade)-2))]
feat.sel.sum$sum <- rowSums(feat.sel.sum)
feat.sel.sum$prevalence <-feat.sel.prev.bacillotaSporeClade$prevalence
feat.sel.sum$present_nPops <-feat.sel.prev.bacillotaSporeClade$present_nPops
# And plot
if (plot_all == TRUE) {
  pdf(file = paste0("/Users/hl636/Desktop/", "sporeCounts_popSpec_", feat, "_", feat_type, "_Bacillota_o_RF39_f_UBA660_final.pdf"), width = 3.3, height = 3)
  print(ggscatter(feat.sel.sum, x = "sum", y = "prevalence", add = "reg.line", size = 1, color = "black", fill="grey60", shape = 21, add.params = list(color = "red2", linetype= "dashed")) + theme_bw() +
          stat_cor(label.x = 5, label.y = 0.2) + xlim(min(feat.sel.sum$sum), max(feat.sel.sum$sum)) + ylim(min(feat.sel.sum$prevalence), max(feat.sel.sum$prevalence)) + ggtitle("Bacillota_o_RF39_f_UBA660: sporulation gene counts") + ylab("prevalence") +
          stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), label.x = 5, label.y = 0.1 ) +
          rremove("grid"))
  print(ggscatter(feat.sel.sum, x = "sum", y = "present_nPops", add = "reg.line", size = 1, color = rgb(0.2,0.2,0.2,0.3), fill=rgb(0.2,0.2,0.2,0.3), shape = 21, add.params = list(color = "red2", linetype= "dashed")) + theme_bw() +
          stat_cor(label.x = min(feat.sel.sum$sum), label.y = min(feat.sel.sum$present_nPops) + 2) + ylim(1,8) + ggtitle("Bacillota_o_RF39_f_UBA660: sporulation gene counts") + ylab("present_nPops") +
          stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), label.x = min(feat.sel.sum$sum), label.y = min(feat.sel.sum$present_nPops) + 1) +
          rremove("grid"))
  dev.off()
}

# Check clade taxonomy
# MAG.metadata3 <- MAG.metadata2
# MAG.metadata3$Order <- str_split_fixed(MAG.metadata3$classification, ";", n = Inf)[,4]
# MAG.metadata3$Family <- str_split_fixed(MAG.metadata3$classification, ";", n = Inf)[,5]
# MAG.metadata3$Order <- gsub("o__", "", MAG.metadata3$Order)
# MAG.metadata3$Family <- gsub("f__", "", MAG.metadata3$Family)
# asd <- merge(bacillota_spore_clade, MAG.metadata3, by = "original_bin")
# qwe <- MAG.metadata3[MAG.metadata3$phylum == "Bacillota",]
# qwe1 <- MAG.metadata3[MAG.metadata3$Order == "RF39",]
# qwe2 <- MAG.metadata3[MAG.metadata3$Family== "UBA660",]
# # Check that all spore-forming clade members belong to the same order/family, and conversely that all members of said order/family are included as spore-formers
# all(qwe2$original_bin == asd$original_bin)