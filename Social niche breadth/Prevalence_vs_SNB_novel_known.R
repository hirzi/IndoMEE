# This script plots prevalence vs SNB

# Load libraries
library(ggplot2)
library(ggpubr)

# Read in data
metadata <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Phylogenetics/HybridAssembly_Comp50Cont5_draft2/Metadata_HybridAssembly_Comp50Cont5_draft2_revised25Nov_MaaslinNeedsToBeUpdated.txt", comment.char="", sep = "\t", header = TRUE)
metadata$novelty <- as.factor(metadata$novelty)

# Plot - basic
print(ggscatter(metadata, x = "prevalence", y = "SNB", add = "reg.line", size = 1, color = rgb(0.2,0.2,0.2,0.2), fill=rgb(0.2,0.2,0.2,0.2), shape = 21, add.params = list(color = "red2", linetype= "dashed")) + theme_bw() +
        stat_cor(label.x = 0.1, label.y = 0.1) + ggtitle("Prevalence vs SNB") + xlim(0,1) + ylim(0,1) +
        stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), label.x = 0.1, label.y = 0.3) +
        rremove("grid"))

# Plot - coloured by novelty
cols <- c("gray50", "lightskyblue")
print(ggscatter(metadata, x = "prevalence", y = "SNB", add = "reg.line", size = 1, color = "novelty", fill="novelty", shape = 21, add.params = list(color = "red2", linetype= "dashed")) + theme_bw() +
        stat_cor(label.x = 0.1, label.y = 0.1) + ggtitle("Prevalence vs SNB") + xlim(0,1) + ylim(0,1) + scale_color_manual(values = cols) + 
        stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), label.x = 0.1, label.y = 0.3) +
        rremove("grid"))
