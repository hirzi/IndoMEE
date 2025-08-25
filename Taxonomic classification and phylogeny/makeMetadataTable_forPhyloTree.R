##### This script produces a metadata table for a specified assembly for use with annotating phylogenetic trees #####

# Load libraries
library(gtools)
library(stringr)
library(randomcoloR)
library(tidyverse)
library(RColorBrewer)

# Set working directories
counts_dir <- "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assemblies/On_hostDecontaminated_fastqs/"
magscreen_dir <- "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/magscreen/Novel_MAGs/"
maaslin_path <- "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/Maaslin/"

# Input assemblies
prefix <- "HybridAssembly_Comp50Cont5_draft2"
magscreen_comp <- "HybridAssembly_sANI95_Comp50Cont5_vs_UHGG_allMAGs"

# List of filtering regimes
filt <- '*total*strong*relAbund*revised.tsv'

# Define input (metamap counts output) files
summary_dir <- paste0(counts_dir, prefix, "/summary")
input <- list.files(summary_dir)[grepl(glob2rx(filt), list.files(summary_dir))]
print(input)

# Read in (relative abudance) data 
counts_table_raw <- read.table(paste0(counts_dir, prefix, "/summary/", input), header = TRUE, sep = "\t")

# Format and filter abundance dataframe
first_sample <- "ASMDAI005"
sample_idx <- match(first_sample, colnames(counts_table_raw))
abund_df <- counts_table_raw[,c(sample_idx:ncol(counts_table_raw))]
non_relevant_samples <- c("POSMOCKCTRL01", "POSMOCKCTRL02")
high_hostDNA_samples <- c("BA012", "BA036", "BAPDW027")
samples_exclude <- c(non_relevant_samples, high_hostDNA_samples)
abund_df <- abund_df[, ! colnames(abund_df) %in% samples_exclude]
n_samples <- ncol(abund_df)

# Calculate meta-population prevalence and abundance
abund_df$totalAbundance <- rowSums(abund_df) / n_samples    
abund_df$Prevalence <- rowSums(abund_df[,c(1:n_samples)]!=0) / n_samples                                        
abund_df <- cbind(counts_table_raw[,c(1,2,3)], abund_df)

# To merge SNB data (which has been deduplicated), we'll need to reformat the counts table toi match
sample_idx <- match(first_sample, colnames(counts_table_raw))
counts_table <- counts_table_raw[, c(3, sample_idx:ncol(counts_table_raw))]
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
# # Merge duplicate entries (add duplicates' counts)
# counts_table_dedup <- aggregate(counts_table[colnames(counts_table)[2:length(colnames(counts_table))]], by=counts_table['classification'], sum)
# counts_table_dedup <- counts_table_dedup[order(counts_table_dedup$classification),]
raw_SNB_df <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Social niche breadth/SNB_matrixOutput_1_.txt", header = TRUE, sep = "\t")

# Get novel MAGs
novel_MAGs <- read.table(paste0(magscreen_dir, magscreen_comp, "/new_species.txt"))$V1
novel_MAGs <- gsub(".fa", "", novel_MAGs)

# Get Maaslin identified population and diet driven MAGs
#pop_MAGs_raw <- read.table(paste0(maaslin_path, "all_results(Sago+Pork+Chicken)_popvar.tsv"), header = TRUE, sep = "\t")
pop_MAGs_raw <- read.table(paste0(maaslin_path, "popMAGs.tsv"), header = TRUE, sep = "\t")
#diet_MAGs_raw <- read.table(paste0(maaslin_path, "all_results(Sago+Pork+Chicken).tsv"), header = TRUE, sep = "\t")
#diet_MAGs_raw <- read.table(paste0(maaslin_path, "significant_results_wTax.tsv"), header = TRUE, sep = "\t")
diet_MAGs_raw <- read.table(paste0(maaslin_path, "dietMAGs.tsv"), header = TRUE, sep = "\t")
diet_MAGs_filt <- diet_MAGs_raw[diet_MAGs_raw$qval < 0.05, ]
diet_MAGs_filt <- diet_MAGs_filt[diet_MAGs_filt$metadata != "Age",]
diet_MAGs_filt <- diet_MAGs_filt[diet_MAGs_filt$metadata != "Gender",]
diet_MAGs_filt <- diet_MAGs_filt[diet_MAGs_filt$metadata != "bmi",]
diet_MAGs <- unique(diet_MAGs_filt[, 1])
pop_MAGs_filt <- pop_MAGs_raw[pop_MAGs_raw$pop_var > 35, ]
pop_MAGs <- pop_MAGs_filt[, 1]
diet_MAGs <- gsub("INDOMEE", "MGYG", diet_MAGs)
pop_MAGs <- gsub("INDOMEE", "MGYG", pop_MAGs)
#pop_MAGs <- sort(counts_table_raw[match(pop_MAGs, counts_table_raw$Genome),][,2])
#diet_MAGs <- sort(counts_table_raw[match(diet_MAGs, counts_table_raw$Genome),][,2])
lifestyle_cline_MAGs_raw <- read.table(paste0(maaslin_path, "LifestyleClineMAGs.tsv"), header = TRUE, sep = "\t")
lifestyle_cline_MAGs_filt <- lifestyle_cline_MAGs_raw[lifestyle_cline_MAGs_raw$qval < 0.05, ]
lifestyle_cline_MAGs_filt <- lifestyle_cline_MAGs_filt[lifestyle_cline_MAGs_filt$metadata != "Age", ]
lifestyle_cline_MAGs_filt <- lifestyle_cline_MAGs_filt[lifestyle_cline_MAGs_filt$metadata != "Gender", ]
lifestyle_cline_MAGs_filt <- lifestyle_cline_MAGs_filt[,c(1,4,9)]
lifestyle_cline_MAGs_filt$feature <- gsub("INDOMEE", "MGYG", lifestyle_cline_MAGs_filt$feature)
colnames(lifestyle_cline_MAGs_filt) <- c("Genome.x", "lifestyle_coef", "lifestyle_qval")

# Get ANI statistics to closest UHGG representative
ANI_df <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/magscreen/Novel_MAGs/HybridAssembly_sANI95_Comp50Cont5_vs_UHGG_allMAGs/summary.tsv", header = TRUE)
colnames(ANI_df)[c(1,2)] <- c("original_bin", "closest_UHGG_ref")

# Get sporulation gene counts from Kegg, GO or pfam
feat <- "kegg"
feat_type <- "binary"
feat.sel.alt <- read.table(paste0("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/Functional enrichment/", feat, "_", feat_type, "_sporeGeneCounts.txt"), sep = "\t", row.names = 1, header = TRUE)
feat.sel.alt$original_bin <- rownames(feat.sel.alt)
feat.sel.alt <- feat.sel.alt[,c(ncol(feat.sel.alt), ncol(feat.sel.alt)-1)]

# Get loadings of MAG abundance PCoA (here calculated with SPMP samples)
pc_loadings1 <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/Functional enrichment/allPhyla_allMAGs_IndowSPMP_MAGAbundPCoA_Aitchison_loadings_revised.txt", header = TRUE, sep = "\t")
pc_loadings1 <- pc_loadings1[,c(1,2)]
pc_loadings1$Genome.x <- rownames(pc_loadings1)
colnames(pc_loadings1) <- c("PC1_Aitchison", "PC2_Aitchison", "Genome.x")
pc_loadings2 <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/Functional enrichment/allPhyla_allMAGs_IndowSPMP_MAGAbundPCoA_Jaccard_loadings_revised.txt", header = TRUE, sep = "\t")
pc_loadings2 <- pc_loadings2[,c(1,2)]
pc_loadings2$Genome.x <- rownames(pc_loadings2)
colnames(pc_loadings2) <- c("PC1_Jaccard", "PC2_Jaccard", "Genome.x")

# Merge metadata
counts_table_metadata <- cbind(counts_table_raw[,c(1,2)], counts_table[,c(1)])
colnames(counts_table_metadata)[ncol(counts_table_metadata)] <- "taxonomic.lineage"
SNB_df <- merge(counts_table_metadata, raw_SNB_df, by = "taxonomic.lineage")
metadata_df <- merge(abund_df, SNB_df, by = "original_bin")
metadata_df <- metadata_df[, c("original_bin", "Genome.x", "classification", "taxonomic.lineage", "totalAbundance", "Prevalence", "number.of.samples", "mean.relative.abundance", "SNB.score")]
metadata_df$Phylum <- str_split_fixed(metadata_df$classification, ";", n = Inf)[,2]
metadata_df$Phylum <- gsub("p__", "", metadata_df$Phylum)
metadata_df$novelty <- as.character(metadata_df$original_bin %in% novel_MAGs)
#metadata_df$novelty <- gsub("TRUE", "Novel", metadata_df$novelty)
metadata_df$novelty <- gsub("TRUE", "1", metadata_df$novelty)
#metadata_df$novelty <- gsub("FALSE", "Known", metadata_df$novelty)
metadata_df$novelty <- gsub("FALSE", "0", metadata_df$novelty)
metadata_df$maaslin_diet <- as.character(metadata_df$Genome.x %in% diet_MAGs)
metadata_df$maaslin_pop <- as.character(metadata_df$Genome.x %in% pop_MAGs)
metadata_df <- metadata_df %>% mutate(maaslin =
                         ifelse(maaslin_diet == TRUE & maaslin_pop == FALSE, "Diet",
                                ifelse(maaslin_diet == FALSE & maaslin_pop == TRUE, "Pop",
                                       ifelse(maaslin_diet == TRUE & maaslin_pop == TRUE, "Diet_Pop",
                                              ifelse(maaslin_diet == FALSE & maaslin_pop == FALSE, NA,
                                                     NA)))))
metadata_df$maaslin_diet <- ifelse(metadata_df$maaslin_diet==TRUE, 1, 0)
metadata_df$maaslin_pop <- ifelse(metadata_df$maaslin_pop==TRUE, 1, 0)
metadata_df <- merge(metadata_df, lifestyle_cline_MAGs_filt, by = "Genome.x", all.x = TRUE)
metadata_df <- merge(metadata_df, ANI_df[,c(1,2,7)], by = "original_bin", all.x = TRUE)
metadata_df <- merge(metadata_df, feat.sel.alt, by = "original_bin", all.x = TRUE)
metadata_df <- merge(metadata_df, pc_loadings1, by = "Genome.x", all.x = TRUE)
metadata_df <- merge(metadata_df, pc_loadings2, by = "Genome.x", all.x = TRUE)
metadata_df[,c(1,2)] <- metadata_df[,c(2,1)]

# Format df
colnames(metadata_df) <- c("original_bin", "genome", "classification", "classification_SNB", "total_abundance", "prevalence", "num_samples", "mean_relative_abundance", "SNB", "phylum", "novelty", "maaslin_diet", "maaslin_pop", "maaslin", "maaslin_lifestyle_coef", "maaslin_lifestyle_qval", "closest_UHGG_ref", "ANI_closest_ref", "n_spore_genes", "PC1_Aitchison", "PC2_Aitchison", "PC1_Jaccard", "PC2_Jaccard")
metadata_df <- metadata_df[, c("original_bin", "genome", "classification", "classification_SNB", "phylum", "closest_UHGG_ref", "ANI_closest_ref", "novelty", "total_abundance", "mean_relative_abundance", "prevalence", "num_samples", "SNB", "maaslin", "maaslin_diet", "maaslin_pop", "maaslin_lifestyle_coef", "maaslin_lifestyle_qval", "n_spore_genes", "PC1_Aitchison", "PC2_Aitchison", "PC1_Jaccard", "PC2_Jaccard")]

# Add colours
metadata_df$phylum <- as.factor(metadata_df$phylum)
#palette <- c("#54FF9F", "#DB5688", "#D75E47", "#E1A863", "#CEC9E1", "#CD0000", "#DAA3D8", "#BBEAD8", "#C99493", "#77C6DF", "#E3DAC2", "#D746E0", "#68228B", "#FFFF00", "#7593D6", "#556B2F", "#6E8B3D", "#A2CD5A", "#CAFF70")
palette <- c("#FFFF66", "#F78375", "#0055AA", "#A65628", "#ACD0F4", "#09979B", "#FFD5E1", "#F781BF", "#CAB2D6", "#984EA3", "#75D8D5", "#E5E5E5", "#26713D", "#9ED470", "#555555", "#FFD700", "#EAD3BF", "#FF7F00", "#C81518") # following Asha's colour palette
metadata_df$phylum <- factor(metadata_df$phylum, levels = c("Thermoplasmatota", "Methanobacteriota", "Bacteroidota", "Verrucomicrobiota", "Elusimicrobiota", "Spirochaetota", "Myxococcota", "Desulfobacterota", "Campylobacterota", "Pseudomonadota", "Synergistota", "Eremiobacterota", "Actinomycetota", "Fusobacteriota", "Cyanobacteriota", "Bacillota", "Bacillota_B", "Bacillota_C", "Bacillota_A"))
#palette <- c("#54FF9F", "#DB5688", "#E3DAC2", "#D746E0", "#68228B", "#DAA3D8", "#BBEAD8", "#C99493", "#77C6DF", "#CEC9E1", "#CD0000", "#E1A863", "#D75E47", "#FFFF00", "#7593D6", "#556B2F", "#6E8B3D", "#A2CD5A", "#CAFF70") # midpoint rooting order (note: this this doesn't change anything here since the phyla will keep it's color via palette[metadata_df$phylum]; however will need to input this order into the iTOL metadata legend description)
#metadata_df$phylum <- factor(metadata_df$phylum, levels = c("Thermoplasmatota", "Methanobacteriota", "Synergistota", "Eremiobacterota", "Actinomycetota", "Myxococcota", "Desulfobacterota", "Campylobacterota", "Pseudomonadota", "Elusimicrobiota", "Spirochaetota", "Verrucomicrobiota", "Bacteroidota", "Fusobacteriota", "Cyanobacteriota", "Bacillota", "Bacillota_B", "Bacillota_C", "Bacillota_A"))  # midpoint rooting order (note: this this doesn't change anything here since the phyla will keep it's color via palette[metadata_df$phylum]; however will need to input this order into the iTOL metadata legend description)
#palette <- distinctColorPalette(length(unique(metadata_df$phylum))) # get random colours
#palette2 <- c("#035887", "#3186B5", "#87CEFF", "#E5E5E5")

# Default rooting
# Midpoint rooting
#levels(metadata_df$phylum)
metadata_df$col_phyum <- palette[metadata_df$phylum]
metadata_df$col_totAbund_idx <- log10(metadata_df$total_abundance)
metadata_df$col_totAbund_idx[is.infinite(metadata_df$col_totAbund_idx)] <- min(metadata_df$col_totAbund_idx[is.finite(metadata_df$col_totAbund_idx)]) - 0.5
totAbund_range <- seq(min(metadata_df$col_totAbund_idx), max(metadata_df$col_totAbund_idx), length.out = 11)
# Use n equally spaced breaks to assign each value to n-1 equal sized bins 
breaks_totAbund <- cut(metadata_df$col_totAbund_idx, breaks = totAbund_range, include.lowest = TRUE)
breaks_prevalence <- cut(metadata_df$prevalence, breaks = seq(0,1,0.1), include.lowest = TRUE)
#breaks_SNB <- cut(metadata_df$SNB, breaks = seq(min(metadata_df$SNB, na.rm = TRUE), max(metadata_df$SNB, na.rm = TRUE), length.out = 11), include.lowest = TRUE)
breaks_SNB <- cut(metadata_df$SNB, breaks = seq(0,1,0.1), include.lowest = TRUE)
breaks_ANI <- cut(metadata_df$ANI_closest_ref, breaks = seq(80,100, length.out = 5), include.lowest = TRUE)
breaks_lifestyleCoef <- cut(metadata_df$maaslin_lifestyle_coef, breaks = seq(-3.3,3.3, length.out = 11), include.lowest = TRUE)
# Use bin indices, ii, to select color from vector of n-1 equally spaced colors
colfunc <- colorRampPalette(c("gray90", "black"))
colfunc2 <- colorRampPalette(c("#87CEFF", "#EEEEEE"))
#metadata_df$col_totAbund <- brewer.pal(10,"RdYlBu")[breaks_totAbund]
#metadata_df$col_prev <- brewer.pal(10,"RdYlBu")[breaks_prevalence]
metadata_df$col_totAbund <- colfunc(10)[breaks_totAbund]
metadata_df$col_prev <- colfunc(10)[breaks_prevalence]
#metadata_df$col_SNB <- brewer.pal(10,"RdYlBu")[breaks_SNB]
metadata_df$col_SNB <- rev(viridis(10, option = "magma"))[breaks_SNB]
metadata_df$col_SNB[is.na(metadata_df$col_SNB)] <- "#FFFFFF"
#metadata_df$col_ANI <- palette2[breaks_ANI]
metadata_df$col_ANI <- colfunc2(4)[breaks_ANI]
metadata_df$col_ANI[is.na(metadata_df$col_ANI)] <- "#FFFFFF"
metadata_df$novel_ANI95 <- ifelse((metadata_df$novelty == "1") & (metadata_df$ANI_closest_ref < 95) & (metadata_df$ANI_closest_ref >= 90), 1, 0)
metadata_df$novel_ANI90 <- ifelse((metadata_df$novelty == "1") & (metadata_df$ANI_closest_ref < 90), 1, 0)
metadata_df$col_lifestyle <- brewer.pal(10,"RdYlBu")[breaks_lifestyleCoef]
metadata_df$col_lifestyle[is.na(metadata_df$col_lifestyle)] <- "#FFFFFF"

#levels(breaks_totAbund)
#brewer.pal(10,"RdYlBu")
# Write out table
write.table(metadata_df, paste0("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Phylogenetics/Metadata_", prefix, "_revised20Aug.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
#feat_df <- unique(metadata_df[,c("phylum", "colour")])
#as.vector(feat_df$phylum)
#as.vector(feat_df$colour)

# Plot x-y scatters
#plot(metadata_df$total_abundance, metadata_df$mean_relative_abundance)
