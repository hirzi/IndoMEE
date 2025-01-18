## This code performs functional enrichment analysis on control-treatment groups (e.g. novel vs non-novel MAGs)

# Load libraries
library(data.table)
library(matrixStats)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(tidyr)
library(ALDEx2)

## Define run parameters
# Run local or remote
run <- "local"
# Choose whether control-treatment groups are known-novel MAGs or population-diet driven MAGs ("novel_known", "pop_diet")
comp <- "novel_known"
# Whether data (features) are counts or binary (i.e. presence-absence) per sample
# binary: gutsmash, antismash, cazy, amrfinder, pfam, GO, Kegg_D
# counts: cog, Kegg_C	
feat_type <- "binary"
# Statistical test. 
if (feat_type == "binary") {
  comp_method <- "fisher"
} else if (feat_type == "counts") {
  # Choose between Wilcoxon rank sum (i.e. Mann-Whitney) test ("wilcoxon"), the same test but via Aldex (which Monte Carlo samples from a Dirichlet distribution to account for uncertainty in counts data, and also calculates effect sizes) ("aldex"), and an Aldex GLM model which includes MAG completeness as a covariate ("aldex_glm")
  comp_method <- "wilcoxon"
}
# Filter to MAGs above certain completeness (only applicable for Wilcoxon and Aldex standard)
compl_filter <- TRUE
compl_theshold <- 90
# Normalise by number of features per MAG (reflective of genome size) - for counts based statistical tests
if (feat_type == "counts") {
  norm_by_feat_length <- TRUE
} else if (feat_type == "binary") {
  norm_by_feat_length <- FALSE
}

# Choose feature dataset (e.g. "cog", "kegg", "gutsmash", "cazy", "antismash", "amrfinder", "pfam", "GO")
feat <- "kegg"
# Define whether to subset to specific phyla
subset_taxa <- FALSE

# Define paths
if (run == "local" | run == "local_remote") {
  feat.path = "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/Functional enrichment/HybridAssembly_Comp50Cont5_Summary/"
  descr.path <- "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/Functional enrichment/"
  maaslin.path <- "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/Maaslin/"
  magscreen.path <- "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/magscreen/Novel_MAGs/"
  all.MAGs <- scan("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/Functional enrichment/dereplicated_genomes.txt", what = "")
  MAG.metadata <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Metadata/metawrap_HybridAssembly_Comp50Cont5_draft2_metadata.tsv", header=TRUE)
  outdir <- "/Users/hl636/Desktop/"
} else if (run == "remote" | run == "remote_local") {
  feat.path = "/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/func_annotation/genofan/Results/HybridAssembly_Comp50Cont5_Summary/"
  descr.path <- "/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/func_analyses/"
  magscreen.path <- "/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/magscreen/Results/"
  all.MAGs <- scan("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_ASSEMBLY/5.4.DREP/Comp50Cont5/dereplicated_genomes.txt", what = "")
  MAG.metadata <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metamap/reference_assemblies_metadata/metawrap_HybridAssembly_Comp50Cont5_draft2_metadata.tsv", header=TRUE)
  outdir <- "/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/func_analyses/"
}

# Load data - features and feature annotations/description
if (feat == "kegg") {
  file <- paste0(feat.path, "kegg_orthologs_summary.tsv")
  no_col <- max(count.fields(file, sep = "\t"))
  feat.data <- read.table(file, sep="\t", fill=TRUE, header = F, col.names=c(1:no_col))
  feat.descr <- read.table(paste0(descr.path, "kegg_descriptions_revised.txt"), header = TRUE, sep = "\t", quote = "")
} else if (feat == "cog") {
  file <- paste0(feat.path, "eggnog_emapper_summary.tsv")
  feat.data = read.delim(file, header=FALSE, stringsAsFactors = FALSE)[,c(1,7)][-1, ]
  feat.descr <- read.table(paste0(descr.path, "cog_categories.txt"), sep = "\t")
  colnames(feat.descr) <- c("feature", "description")
} else if (feat == "gutsmash") {
  file <- paste0(feat.path, "gutsmash_results_summary.tsv")
  feat.data <- read.table(file, header=TRUE)
  feat.data <- feat.data[,c(1,3)]
} else if (feat == "cazy") {
  file <- paste0(feat.path, "cazy_results_summary.tsv")
  no_col <- max(count.fields(file, sep = "\t"))
  feat.data <- read.table(file, sep="\t", fill=TRUE, header = F, col.names=c(1:no_col))
  #feat.descr <- read.table(paste0(descr.path, ""), header = TRUE, sep = "\t", quote = "") # need to look for list of cazy terms
} else if (feat == "antismash") {
  file <- paste0(feat.path, "antismash_results_summary.tsv")
  feat.data <- read.table(file, header=TRUE)
  feat.data$genome <- gsub("gnlX", "", feat.data$genome)
} else if (feat == "amrfinder") {
  file <- paste0(feat.path, "amrfinder_results_summary.tsv")
  feat.data <- read.table(file, header = TRUE, sep = "\t", quote = "")
  feat.data$Contig.id <- gsub("gnl\\|X\\|", "", feat.data$Contig.id)
  feat.data <- feat.data[,c("Contig.id", "Gene.symbol")] # Choose whether to run tests on gene or sub/class level
  #feat.data <- feat.data[,c("Contig.id", "Subclass")]
  #feat.data <- feat.data[,c("Contig.id", "Class")]
} else if (feat == "GO") {
  file <- paste0(feat.path, "eggnog_emapper_summary.tsv")
  feat.data = read.delim(file, header=FALSE, stringsAsFactors = FALSE)[,c(1,10)][-1, ]
  feat.descr <- read.table(paste0(descr.path, "GO_descriptions.txt"), header = TRUE, sep = "\t", comment.char = "", quote = "")
} else if (feat == "pfam") {
  file <- paste0(feat.path, "eggnog_emapper_summary.tsv")
  feat.data = read.delim(file, header=FALSE, stringsAsFactors = FALSE)[,c(1,21)][-1, ]
}

# Load data - MAGs
if (comp == "novel_known") {
  magscreen_comp <- "HybridAssembly_sANI95_Comp50Cont5_vs_UHGG_allMAGs"
  novel.MAGs <- scan(paste0(magscreen.path, magscreen_comp, "/new_species.txt"), what = "")
  novel.MAGs <- gsub(".fa", "", novel.MAGs)
  all.MAGs <- gsub(".fa", "", all.MAGs)
  known.MAGs <- all.MAGs[which(!all.MAGs %in% novel.MAGs)]
  #known.MAGs <- sample(known.MAGs, length(novel.MAGs), replace = FALSE)
  control.MAGs <- known.MAGs
  treatment.MAGs <- novel.MAGs
} else if (comp == "pop_diet") {
  #pop_MAGs_raw <- read.table(paste0(maaslin.path, "best_population_effects.tsv"), header = TRUE, sep = "\t")
  pop_MAGs_raw <- read.table(paste0(maaslin.path, "all_results(Sago+Pork+Chicken)_popvar.tsv"), header = TRUE, sep = "\t")
  #diet_MAGs_raw <- read.table(paste0(maaslin.path, "diet-only_MAGs.tsv"), header = TRUE, sep = "\t")
  #diet_MAGs_raw <- read.table(paste0(maaslin.path, "all_results(Sago+Pork+Chicken).tsv"), header = TRUE, sep = "\t")
  diet_MAGs_raw <- read.table(paste0(maaslin.path, "significant_results_wTax.tsv"), header = TRUE, sep = "\t")
  pop_MAGs_raw <- pop_MAGs_raw[order(pop_MAGs_raw$pop_var, decreasing = TRUE),]  
  diet_MAGs_raw <- diet_MAGs_raw[order(diet_MAGs_raw$qval),]
  diet_MAGs_filt <- diet_MAGs_raw[diet_MAGs_raw$qval<0.05,]
  pop_MAGs_filt <- pop_MAGs_raw[pop_MAGs_raw$pop_var>35, ]
  diet_MAGs <- unique(diet_MAGs_filt[, 1])
  pop_MAGs <- pop_MAGs_filt$feature
  if (length(pop_MAGs) > length(diet_MAGs)) {
    pop_MAGs <- pop_MAGs[1:length(diet_MAGs)]
  } else if (length(pop_MAGs) < length(diet_MAGs)) {
    diet_MAGs <- diet_MAGs[1:length(pop_MAGs)]
  }
  diet_pop_MAGs_intersect <- intersect(pop_MAGs, diet_MAGs) # remove intersecting elements, to make sets mutually exclusive
  diet_MAGs <- diet_MAGs [! diet_MAGs %in% diet_pop_MAGs_intersect]
  pop_MAGs <- pop_MAGs [! pop_MAGs %in% diet_pop_MAGs_intersect]
  #control.MAGs <- sort(MAG.metadata[match(pop_MAGs, MAG.metadata$Genome),][,2])
  #treatment.MAGs <- sort(MAG.metadata[match(diet_MAGs, MAG.metadata$Genome),][,2])
  control.MAGs <- pop_MAGs
  treatment.MAGs <- diet_MAGs
}

# Subset data to specific phyla if applicable
if (subset_taxa == TRUE) {
  sel.taxa <- paste0("p__", c("Actinomycetota", "Bacillota", "Pseudomonadota", "Mycoplasmatota"))
  tax.annot.df <- as.data.frame(str_split_fixed(MAG.metadata$classification, ";", 7))
  colnames(tax.annot.df) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  MAG.metadata <- cbind(MAG.metadata, tax.annot.df)
  MAG.metadata <- MAG.metadata[MAG.metadata$phylum %in% sel.taxa, ]
  all.MAGs <- all.MAGs[which(all.MAGs %in% MAG.metadata$original_bin)]
  treatment.MAGs <- treatment.MAGs[which(treatment.MAGs %in% MAG.metadata$original_bin)]
  control.MAGs <- control.MAGs[which(control.MAGs %in% MAG.metadata$original_bin)]
}

# Filter MAGs if applicable
if (compl_filter == TRUE) {
  if (comp_method == "wilcoxon" | comp_method == "aldex") {
    MAG.metadata.filt <- MAG.metadata[MAG.metadata$completeness > compl_theshold, ]
    all.MAGs <- all.MAGs[which(all.MAGs %in% MAG.metadata.filt$original_bin)]
    treatment.MAGs <- treatment.MAGs[which(treatment.MAGs %in% MAG.metadata.filt$original_bin)]
    control.MAGs <- control.MAGs[which(control.MAGs %in% MAG.metadata.filt$original_bin)]
  }
}

# Extract features
if (feat == "kegg" | feat == "cog" | feat == "cazy" | feat == "amrfinder" | feat == "GO" | feat == "pfam") {
  proteins.prefixes <-  sub("_[^_]*$", "", feat.data[,1])
} else if (feat == "gutsmash" | feat == "antismash") {
  proteins.prefixes <-  feat.data[,1]
}
i <- 1
# For treatment MAGs
for (MAG in treatment.MAGs) {
  feat.treatment.MAG = as.vector(as.matrix(feat.data[which(proteins.prefixes %in% MAG),-1]))
  feat.treatment.MAG = feat.treatment.MAG[which(feat.treatment.MAG != "")]
  if (feat == "kegg") {
    # Choose below between function_B, function_C and function_D
    # Note that while each KEGG has a single function_D, it can have multiple function_C and function_B's (i.e. the same gene can be involved in multiple pathways). Here, we allow single genes' contribution to multiple pathways (i.e. # of features can be > than number of keggs)
    feat.treatment.MAG <- feat.descr[which(feat.descr$kegg %in% feat.treatment.MAG),]$function_D
  } else if (feat == "cog") {
    multi_feats <- feat.treatment.MAG[nchar(feat.treatment.MAG) > 1]
    single_feats <- feat.treatment.MAG[nchar(feat.treatment.MAG) == 1]
    multi_feat.unlist <- unlist(sapply(multi_feats, function(x) unlist(strsplit(x, split = ""))))
    feat.treatment.MAG <- c(single_feats, multi_feat.unlist)
  } else if (feat == "GO" | feat == "pfam") {
    feat.treatment.MAG = feat.treatment.MAG[which(feat.treatment.MAG != "-")]
    feat.treatment.MAG <- unlist(strsplit(feat.treatment.MAG, ","))
    if (feat == "GO") {
      feat.treatment.MAG <- gsub("GO:", "", feat.treatment.MAG)
    }
  }
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
  if (feat == "kegg") {
    # Choose below between function_B, function_C and function_D
    # Note that while each KEGG has a single function_D, it can have multiple function_C and function_B's (i.e. the same gene can be involved in multiple pathways). Here, we allow single genes' contribution to multiple pathways (i.e. # of features can be > than number of keggs)
    feat.control.MAG <- feat.descr[which(feat.descr$kegg %in% feat.control.MAG),]$function_D
  } else if (feat == "cog") {
    multi_feats <- feat.control.MAG[nchar(feat.control.MAG) > 1]
    single_feats <- feat.control.MAG[nchar(feat.control.MAG) == 1]
    multi_feat.unlist <- unlist(sapply(multi_feats, function(x) unlist(strsplit(x, split = ""))))
    feat.control.MAG <- c(single_feats, multi_feat.unlist)
  } else if (feat == "GO" | feat == "pfam") {
    feat.control.MAG = feat.control.MAG[which(feat.control.MAG != "-")]
    feat.control.MAG <- unlist(strsplit(feat.control.MAG, ","))
    if (feat == "GO") {
      feat.control.MAG <- gsub("GO:", "", feat.control.MAG)
    }
  }
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

# Normalise by number of features
if (norm_by_feat_length == TRUE) {
  feat.all <- apply(feat.all, 2, function(x) {x/sum(x)})
  feat.all[is.na(feat.all)] <- 0
}
    
# Aldex requires input to be integers, so transform accordingly
if (comp_method == "aldex" || comp_method == "aldex_glm") {
  feat.all <- round(feat.all * 1/min(feat.all[feat.all > 0]))
}
feat.all.original <- feat.all

# Test for significance of difference (in feature counts) between the two groups (treatment vs control MAGs)
if (comp_method == "wilcoxon") {
  # Run Wilcoxon test to test for difference between 2 groups
  k.pvalues <- c()
  k.names <- c()
  k.means.treatment <- c()
  k.means.control <- c()
  for (k in seq(1, nrow(feat.all))) {
    treatment.k.vec <- as.vector(t(feat.all[k,treatment_idx]))
    control.k.vec <- as.vector(t(feat.all[k,control_idx]))
    k.pvalue <- wilcox.test(treatment.k.vec, control.k.vec, paired = FALSE)$p.value
    k.pvalues <- c(k.pvalues, k.pvalue)
    k.names <- c(k.names, row.names(feat.all)[k])
    k.means.treatment <- c(k.means.treatment, mean(treatment.k.vec))
    k.means.control <- c(k.means.control, mean(control.k.vec))
  }
  k.fdr <- p.adjust(k.pvalues, method="fdr")
  res.all <- data.frame(k.names, k.pvalues, k.fdr, k.means.treatment, k.means.control)
  colnames(res.all) <- c("feature", "p_value", "wi.eBH", "mean_treatment", "mean_control")
  res.all$mean_diff <- (res.all$mean_treatment - res.all$mean_control) / res.all$mean_treatment
  
} else if (comp_method == "aldex" || comp_method == "aldex_glm") {
  # Use Aldex to test for difference between 2 groups. For reference, see: https://www.bioconductor.org/packages/release/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.html
  # Aldex can be slow for large datasets, so run on cluster
  if (run == "remote" | run == "local") {
    k.groups = as.factor(c(rep("treatment",length(treatment_idx)),rep("control", length(control_idx))))
    if (comp_method == "aldex") {
      aldex.analy = aldex.clr(reads=feat.all, conds=k.groups, mc.samples=128, denom="all", verbose=TRUE, useMC=TRUE) # ideally go for 128 but would require running on cluster due to memory limits (Error: vector memory exhausted (limit reached?))
      aldex.eff = aldex.effect(aldex.analy, useMC=TRUE, CI=TRUE, verbose=TRUE) #identify those features that where the 95% CI of the effect size does not cross 0.
      aldex.res = aldex.ttest(aldex.analy, verbose=TRUE)
      res.all = data.frame(rownames(aldex.eff), aldex.eff,aldex.res)
    } else if (comp_method == "aldex_glm") {
      # Alternatively, add genome completeness as a covariate. For reference, see: https://genomicsclass.github.io/book/pages/expressing_design_formula.html
      #k.groups <- relevel(k.groups, "treatment") # switch reference level
      MAG.metadata <- MAG.metadata[match(gsub("^.{0,6}", "", colnames(feat.all)), MAG.metadata$original_bin),] # reorder metadata to match order of colnames(keg.all)
      covariates <- data.frame("group" = k.groups, "completeness" =  MAG.metadata$completeness)
      mm <- model.matrix(~ group + completeness, covariates)
      aldex.analy <- aldex.clr(reads=feat.all, mm, mc.samples=128, denom="all", verbose=TRUE, useMC=TRUE)
      aldex.res <- aldex.glm(aldex.analy, mm)
      aldex.eff <- aldex.glm.effect(aldex.analy, useMC=TRUE, CI=TRUE, verbose=TRUE)
      res.all = data.frame(rownames(aldex.eff$group), aldex.eff$group, aldex.res)
    }
  } else if (run == "local_remote") {
    aldex.analy <- readRDS(paste0(outdir, comp_method, "_analy_filter", as.character(compl_filter), compl_theshold, ".rds"))
    aldex.eff <- readRDS(paste0(outdir, comp_method, "_eff_filter", as.character(compl_filter), compl_theshold, ".rds"))
    aldex.res <- readRDS(paste0(outdir, comp_method, "_res_filter", as.character(compl_filter), compl_theshold, ".rds"))
    res.all <- readRDS(paste0(outdir, comp_method, "_res_all_filter", as.character(compl_filter), compl_theshold, ".rds"))
  }
} else if (comp_method == "fisher") {
  # Convert counts to binary presence-absence
  feat.all[feat.all > 1] <- 1
  treatment_absent <- rowSums(feat.all[,treatment_idx] == 0)
  treatment_present <- rowSums(feat.all[,treatment_idx] == 1)
  control_absent <- rowSums(feat.all[,control_idx] == 0)
  control_present <- rowSums(feat.all[,control_idx] == 1)
  cont_df <- cbind(treatment_present, treatment_absent, control_present, control_absent)
  cont_df <- round(cont_df * (1 / min(cont_df[cont_df > 0])))
  k.pvalues <- c()
  k.names <- c()
  k.odds_ratios <- c()
  for (k in seq(1, nrow(cont_df))) {
    cont_feat <- cont_df[k,]
    cont_table <- matrix(c(cont_feat[1], cont_feat[3], cont_feat[2], cont_feat[4]), nrow = 2,
                          dimnames = list(c("Treatment", "Control"), c("Presence", "Absence")))
    k.names <- c(k.names, row.names(cont_df)[k])
    k.fisher <- fisher.test(cont_table)
    k.pvalue <- k.fisher$p.value
    k.odds_ratio <- k.fisher$estimate
    k.pvalues <- c(k.pvalues, k.pvalue)
    k.odds_ratios <- c(k.odds_ratios, k.odds_ratio)
  }
  k.fdr <- p.adjust(k.pvalues, method="fdr")
  res.all <- data.frame(k.names, k.pvalues, k.fdr, k.odds_ratios)
  colnames(res.all) <- c("feature", "p_value", "wi.eBH", "odds_ratio")
}
    
# Save/read in aldex output (since aldex is resource intensive and takes a long time to run, run it on the cluster)
if (run == "remote_local") {
  if (comp_method == "aldex" | comp_method == "aldex_glm") {
    saveRDS(aldex.analy, paste0(outdir, comp_method, "_analy_filter", as.character(compl_filter), compl_theshold, ".rds"))
    saveRDS(aldex.eff, paste0(outdir, comp_method, "_eff_filter", as.character(compl_filter), compl_theshold, ".rds"))
    saveRDS(aldex.res, paste0(outdir, comp_method, "_res_filter", as.character(compl_filter), compl_theshold, ".rds"))
    saveRDS(res.all, paste0(outdir, comp_method, "_res_all_filter", as.character(compl_filter), compl_theshold, ".rds"))
  }
}

if (run == "local" | run == "local_remote") {
  
  # Check if any features DONT overlap 0
  if (comp_method == "aldex") {
    res.all$overlap_zero = ifelse(res.all$effect.low <= 0 & 0 <= res.all$effect.high, "TRUE", "FALSE")
  }
  
  # Summarize results and get most significant
  colnames(res.all)[1] <- "feature"
  if (comp_method == "wilcoxon" || comp_method == "aldex" || comp_method == "fisher") {
    res.sig <- res.all[which(res.all$wi.eBH < 0.05),]
  } else if (comp_method == "aldex_glm") {
    res.sig <- res.all[which(res.all$grouptreatment.pval.holm < 0.05),]
  }
  if(nrow(res.sig) > 0) {
    
    if (comp_method == "aldex" || comp_method == "aldex_glm") {
      effect_thresh = 0.1
      res.sig <- res.sig[which(abs(res.sig$effect) > effect_thresh),]
    } else if (comp_method == "wilcoxon") {
      table(res.sig$direction)
    }
    if (feat == "cog") {
      res.sig <- merge(res.sig, feat.descr, by="feature")
      res.plot <- res.sig
    } else if (feat == "GO") {
      feat.descr$id <- gsub("GO:", "", feat.descr$id)
      res.sig <- merge(res.sig, feat.descr, by.x = "feature", by.y = "id")
      res.sig$description <- paste0(res.sig$feature, ": ", res.sig$name)
    }
    
    # Plot barcharts with effect size
    colnames(res.sig)[match("wi.eBH", colnames(res.sig))] <- "p.adjust"
    colnames(res.sig)[match("odds_ratio", colnames(res.sig))] <- "effect_size"
    colnames(res.sig)[match("mean_diff", colnames(res.sig))] <- "effect_size"
    res.sig = res.sig[order(res.sig$effect_size, decreasing=TRUE),]
    if (feat == "kegg") {
      res.sig$description <- str_sub(res.sig$feature, 6, -15)
    } else if (feat == "pfam") {
      res.sig$description <- res.sig$feature
    }
    if (comp_method == "fisher"){
      if (feat != "GO") {
        res.sig$description <- res.sig$feature
      }
      res.sig$effect_size <- log10(res.sig$effect_size)
      y.lab = "log10(Odds Ratio)"
    } else {
      y.lab = "Effect size"
    }
    # Deal with infinite effect size values
    inf_neg_idx <- which(is.infinite(res.sig$effect_size) & res.sig$effect_size < 0)
    inf_pos_idx <- which(is.infinite(res.sig$effect_size) & res.sig$effect_size > 0)
    if (length(inf_neg_idx) > 0 | length(inf_pos_idx) > 0) {
      #inf_neg_val <- floor(min(res.sig$effect_size[is.finite(res.sig$effect_size)]))
      inf_neg_val <- min(res.sig$effect_size[is.finite(res.sig$effect_size)])
      #inf_pos_val <- ceiling(max(res.sig$effect_size[is.finite(res.sig$effect_size)]))
      inf_pos_val <- max(res.sig$effect_size[is.finite(res.sig$effect_size)])
      res.sig$effect_size[res.sig$effect_size == -Inf] <- inf_neg_val
      res.sig$effect_size[res.sig$effect_size == Inf] <- inf_pos_val
      inf_idx <- c(inf_neg_idx, inf_pos_idx)
      inf_vals <- c(rep(inf_neg_val,length(inf_neg_idx)), rep(inf_pos_val,length(inf_pos_idx)))
      inf_df <- data.frame(inf_idx, inf_vals)
    } else {
      inf_df <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("inf_idx", "inf_vals"))
    }
    # Plot all significant features
    print(ggplot()
          + geom_bar(data=res.sig, aes(x=description, y=effect_size, fill=p.adjust), stat="identity", alpha=0.8, size=0.4)
          + geom_point(data=inf_df, aes(y = inf_vals, x = inf_idx), position = position_dodge(width = 0.9), shape=8, size=5, show.legend=FALSE)
          + theme_bw()
          + rremove("grid")
          + ggtitle(paste0(feat, "_", feat_type))
          + coord_flip()
          + theme(plot.title = element_text(hjust = 0.5))
          + ylab(y.lab)
          + scale_x_discrete(limits=res.sig$description)
          + scale_fill_gradient(low = "indianred1", high = "red3", na.value = NA, trans = 'reverse')
          + theme(axis.title.y = element_blank()))
    
    # Plot top N signficant features, plus features of interest
    top_n <- 10
    topn.pos <- res.sig[res.sig$effect_size > 0, ][1:top_n,]
    topn.neg <- tail(res.sig[res.sig$effect_size < 0, ], top_n)
    if (feat == "GO" | feat == "kegg") {
      res.sig.spore1 <- res.sig[grep("spore", res.sig$description), ]
      res.sig.spore2 <- res.sig[grep("sporu", res.sig$description), ]
    } else if (feat == "pfam") {
      res.sig.spore1 <- res.sig[grep("Spo", res.sig$description), ]
    }
    res.sig.top <- merge(topn.pos, topn.neg, all=TRUE)
    res.sig.top <- merge(res.sig.top, res.sig.spore1, all=TRUE)
    if (feat == "GO" | feat == "kegg") {
      res.sig.top <- merge(res.sig.top, res.sig.spore2, all=TRUE)
    }
    res.sig.top <- na.omit(res.sig.top)
    res.sig.top <- res.sig.top[order(res.sig.top$effect_size, decreasing=TRUE),]
    
    inf_neg_idx2 <- which(is.infinite(res.sig.top$effect_size) & res.sig.top$effect_size < 0)
    inf_pos_idx2 <- which(is.infinite(res.sig.top$effect_size) & res.sig.top$effect_size > 0)
    if (length(inf_neg_idx2) > 0 | length(inf_pos_idx2) > 0) {
      inf_neg_val2 <- floor(min(res.sig.top$effect_size[is.finite(res.sig.top$effect_size)]))
      inf_pos_val2 <- ceiling(max(res.sig.top$effect_size[is.finite(res.sig.top$effect_size)]))
      res.sig.top$effect_size[res.sig.top$effect_size == -Inf] <- inf_neg_val2
      res.sig.top$effect_size[res.sig.top$effect_size == Inf] <- inf_pos_val2
      inf_idx2 <- c(inf_neg_idx2, inf_pos_idx2)
      inf_vals2 <- c(rep(inf_neg_val2,length(inf_neg_idx2)), rep(inf_pos_val2,length(inf_pos_idx2)))
      inf_df_topn <- data.frame(inf_idx2, inf_vals2)
      colnames(inf_df_topn) <- c("inf_idx", "inf_vals")
    } else {
      inf_df_topn <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("inf_idx", "inf_vals"))
    }
    # Add short descriptor
    res.sig.top$description_short <- str_split_fixed(res.sig.top$description, ";", n = Inf)[,1]
    # Plot
    print(ggplot()
          + geom_bar(data=res.sig.top, aes(x=description_short, y=effect_size, fill=p.adjust), stat="identity", alpha=0.8, size=0.4)
          + geom_point(data=inf_df_topn, aes(y = inf_vals, x = inf_idx), position = position_dodge(width = 0.9), shape=8, size=5, show.legend=FALSE)
          + theme_bw()
          + rremove("grid")
          + ggtitle(paste0(feat, "_", feat_type))
          + coord_flip()
          + theme(plot.title = element_text(hjust = 0.5))
          + ylab(y.lab)
          + scale_x_discrete(limits=res.sig.top$description_short)
          + scale_fill_gradient(low = "indianred1", high = "red3", na.value = NA, trans = 'reverse')
          + theme(axis.title.y = element_blank()))
  } else {
    print("No signicant results")
  }
}

# Write out results table
#write.table(res.sig, paste0("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/Functional enrichment/", feat, "_", feat_type, ".txt"), quote = FALSE, row.names = FALSE, sep = "\t")



####################################################################################################

# Plot top N significant features (only)
top_n <- 10
topn.pos <- res.sig[res.sig$effect_size > 0, ][1:top_n,]
topn.neg <- tail(res.sig[res.sig$effect_size < 0, ], top_n)
res.sig.top <- merge(topn.pos, topn.neg, all=TRUE)
res.sig.top <- na.omit(res.sig.top)
res.sig.top <- res.sig.top[order(res.sig.top$effect_size, decreasing=TRUE),]
inf_neg_idx2 <- which(is.infinite(res.sig.top$effect_size) & res.sig.top$effect_size < 0)
inf_pos_idx2 <- which(is.infinite(res.sig.top$effect_size) & res.sig.top$effect_size > 0)
if (length(inf_neg_idx2) > 0 | length(inf_pos_idx2) > 0) {
  inf_neg_val2 <- floor(min(res.sig.top$effect_size[is.finite(res.sig.top$effect_size)]))
  inf_pos_val2 <- ceiling(max(res.sig.top$effect_size[is.finite(res.sig.top$effect_size)]))
  res.sig.top$effect_size[res.sig.top$effect_size == -Inf] <- inf_neg_val2
  res.sig.top$effect_size[res.sig.top$effect_size == Inf] <- inf_pos_val2
  inf_idx2 <- c(inf_neg_idx2, inf_pos_idx2)
  inf_vals2 <- c(rep(inf_neg_val2,length(inf_neg_idx2)), rep(inf_pos_val2,length(inf_pos_idx2)))
  inf_df_topn <- data.frame(inf_idx2, inf_vals2)
  colnames(inf_df_topn) <- c("inf_idx", "inf_vals")
} else {
  inf_df_topn <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("inf_idx", "inf_vals"))
}
# Add short descriptor
res.sig.top$description_short <- str_split_fixed(res.sig.top$description, ";", n = Inf)[,1]
# Annotate sporing MAGs
if (feat == "GO" | feat == "kegg") {
  idx_sporing <- sort(c(grep("spore", res.sig.top$description), grep("sporu", res.sig.top$description)))
} else if (feat == "pfam") {
  idx_sporing <- grep("Spo", res.sig.top$description)
}
res.sig.top[, ncol(res.sig.top) + 1] <- 0
colnames(res.sig.top)[ncol(res.sig.top)] <- "sporing"
res.sig.top$sporing[idx_sporing] <- 1
res.sig.top$sporing <- as.factor(res.sig.top$sporing)

# Plot
print(ggplot()
      + geom_bar(data=res.sig.top, aes(x=description_short, y=effect_size, fill=sporing), stat="identity", alpha=0.8, size=0.4)
      + geom_point(data=inf_df_topn, aes(y = inf_vals, x = inf_idx), position = position_dodge(width = 0.9), shape=8, size=5, show.legend=FALSE)
      + theme_bw()
      + rremove("grid")
      + ggtitle(paste0(feat, "_", feat_type))
      + coord_flip()
      + theme(plot.title = element_text(hjust = 0.5))
      + ylab(y.lab)
      + scale_x_discrete(limits=res.sig.top$description_short)
      + scale_fill_manual(values = c("red3", "dodgerblue3"))
      + theme(axis.title.y = element_blank()))

###############################################################

# Plot features of interest (i.e. sporing MAGs) only
if (feat == "GO" | feat == "kegg") {
  res.sig.spore <- res.sig[grep("spore", res.sig$description), ]
  res.sig.spore2 <- res.sig[grep("sporu", res.sig$description), ]
} else if (feat == "pfam") {
  res.sig.spore <- res.sig[grep("Spo", res.sig$description), ]
}
if (feat == "GO" | feat == "kegg") {
  res.sig.spore <- merge(res.sig.spore, res.sig.spore2, all=TRUE)
}
res.sig.spore <- na.omit(res.sig.spore)
res.sig.spore <- res.sig.spore[order(res.sig.spore$effect_size, decreasing=TRUE),]
inf_neg_idx2 <- which(is.infinite(res.sig.spore$effect_size) & res.sig.spore$effect_size < 0)
inf_pos_idx2 <- which(is.infinite(res.sig.spore$effect_size) & res.sig.spore$effect_size > 0)
if (length(inf_neg_idx2) > 0 | length(inf_pos_idx2) > 0) {
  inf_neg_val2 <- floor(min(res.sig.spore$effect_size[is.finite(res.sig.spore$effect_size)]))
  inf_pos_val2 <- ceiling(max(res.sig.spore$effect_size[is.finite(res.sig.spore$effect_size)]))
  res.sig.spore$effect_size[res.sig.spore$effect_size == -Inf] <- inf_neg_val2
  res.sig.spore$effect_size[res.sig.spore$effect_size == Inf] <- inf_pos_val2
  inf_idx2 <- c(inf_neg_idx2, inf_pos_idx2)
  inf_vals2 <- c(rep(inf_neg_val2,length(inf_neg_idx2)), rep(inf_pos_val2,length(inf_pos_idx2)))
  inf_df_topn <- data.frame(inf_idx2, inf_vals2)
  colnames(inf_df_topn) <- c("inf_idx", "inf_vals")
} else {
  inf_df_topn <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("inf_idx", "inf_vals"))
}
# Add short descriptor
res.sig.spore$description_short <- str_split_fixed(res.sig.spore$description, ";", n = Inf)[,1]
res.sig.spore[, ncol(res.sig.spore) + 1] <- 1
colnames(res.sig.spore)[ncol(res.sig.spore)] <- "sporing"
res.sig.spore$sporing <- as.factor(res.sig.spore$sporing)

# Plot
print(ggplot()
      + geom_bar(data=res.sig.spore, aes(x=description_short, y=effect_size, fill=sporing), stat="identity", alpha=0.8, size=0.4)
      + geom_point(data=inf_df_topn, aes(y = inf_vals, x = inf_idx), position = position_dodge(width = 0.9), shape=8, size=5, show.legend=FALSE)
      + theme_bw()
      + rremove("grid")
      + ggtitle(paste0(feat, "_", feat_type))
      + coord_flip()
      + theme(plot.title = element_text(hjust = 0.5))
      + ylab(y.lab)
      + scale_x_discrete(limits=res.sig.spore$description_short)
      #+ scale_fill_gradient(low = "indianred1", high = "red3", na.value = NA, trans = 'reverse')
      + scale_fill_manual(values = c("dodgerblue3"))
      + theme(axis.title.y = element_blank()))

###############################################################

# Plot top N signficant features, plus features of interest; colouring by feature
top_n <- 10
topn.pos <- res.sig[res.sig$effect_size > 0, ][1:top_n,]
topn.neg <- tail(res.sig[res.sig$effect_size < 0, ], top_n)
if (feat == "GO" | feat == "kegg") {
  res.sig.spore1 <- res.sig[grep("spore", res.sig$description), ]
  res.sig.spore2 <- res.sig[grep("sporu", res.sig$description), ]
} else if (feat == "pfam") {
  res.sig.spore1 <- res.sig[grep("Spo", res.sig$description), ]
}
res.sig.top <- merge(topn.pos, topn.neg, all=TRUE)
res.sig.top <- merge(res.sig.top, res.sig.spore1, all=TRUE)
if (feat == "GO" | feat == "kegg") {
  res.sig.top <- merge(res.sig.top, res.sig.spore2, all=TRUE)
}
res.sig.top <- na.omit(res.sig.top)
res.sig.top <- res.sig.top[order(res.sig.top$effect_size, decreasing=TRUE),]

inf_neg_idx2 <- which(is.infinite(res.sig.top$effect_size) & res.sig.top$effect_size < 0)
inf_pos_idx2 <- which(is.infinite(res.sig.top$effect_size) & res.sig.top$effect_size > 0)
if (length(inf_neg_idx2) > 0 | length(inf_pos_idx2) > 0) {
  inf_neg_val2 <- floor(min(res.sig.top$effect_size[is.finite(res.sig.top$effect_size)]))
  inf_pos_val2 <- ceiling(max(res.sig.top$effect_size[is.finite(res.sig.top$effect_size)]))
  res.sig.top$effect_size[res.sig.top$effect_size == -Inf] <- inf_neg_val2
  res.sig.top$effect_size[res.sig.top$effect_size == Inf] <- inf_pos_val2
  inf_idx2 <- c(inf_neg_idx2, inf_pos_idx2)
  inf_vals2 <- c(rep(inf_neg_val2,length(inf_neg_idx2)), rep(inf_pos_val2,length(inf_pos_idx2)))
  inf_df_topn <- data.frame(inf_idx2, inf_vals2)
  colnames(inf_df_topn) <- c("inf_idx", "inf_vals")
} else {
  inf_df_topn <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("inf_idx", "inf_vals"))
}
# Add short descriptor
res.sig.top$description_short <- str_split_fixed(res.sig.top$description, ";", n = Inf)[,1]
# Annotate sporing MAGs
if (feat == "GO" | feat == "kegg") {
  idx_sporing <- sort(c(grep("spore", res.sig.top$description), grep("sporu", res.sig.top$description)))
} else if (feat == "pfam") {
  idx_sporing <- grep("Spo", res.sig.top$description)
}
res.sig.top[, ncol(res.sig.top) + 1] <- 0
colnames(res.sig.top)[ncol(res.sig.top)] <- "sporing"
res.sig.top$sporing[idx_sporing] <- 1
res.sig.top$sporing <- as.factor(res.sig.top$sporing)

# Plot
print(ggplot()
      + geom_bar(data=res.sig.top, aes(x=description_short, y=effect_size, fill=sporing), stat="identity", alpha=0.8, size=0.4)
      + geom_point(data=inf_df_topn, aes(y = inf_vals, x = inf_idx), position = position_dodge(width = 0.9), shape=8, size=5, show.legend=FALSE)
      + theme_bw()
      + rremove("grid")
      + ggtitle(paste0(feat, "_", feat_type))
      + coord_flip()
      + theme(plot.title = element_text(hjust = 0.5))
      + ylab(y.lab)
      + scale_x_discrete(limits=res.sig.top$description_short)
      + scale_fill_manual(values = c("red3", "dodgerblue3"))
      + theme(axis.title.y = element_blank()))

###############################################################

# Plot all significant features, highlighting sporulation-associated genes
res.sig.plot <- res.sig
# Add short descriptor
res.sig.plot$description_short <- str_split_fixed(res.sig.plot$description, ";", n = Inf)[,1]
# Rename genes incase duplicate names remain after shortening
n_occur <- data.frame(table(res.sig.plot$description_short))
idx_dup_i <- which(n_occur$Freq > 1)
for (i in idx_dup_i) {
  gene_name <- as.character(n_occur$Var1[i])
  gene_dup_freq <- n_occur$Freq[i]
  idx_dup_j <- which(res.sig.plot$description_short == gene_name)
  k <- 1
  for (j in idx_dup_j) {
    res.sig.plot$description_short[j] <- paste0(gene_name, "_", k)
    k <- k + 1
  }
}
# Annotate sporing MAGs
if (feat == "GO" | feat == "kegg") {
  idx_sporing <- sort(c(grep("spore", res.sig.plot$description), grep("sporu", res.sig.plot$description)))
} else if (feat == "pfam") {
  idx_sporing <- grep("Spo", res.sig.plot$description)
}
res.sig.plot[, ncol(res.sig.plot) + 1] <- 0
colnames(res.sig.plot)[ncol(res.sig.plot)] <- "sporing"
res.sig.plot$sporing[idx_sporing] <- 1
res.sig.plot$sporing <- as.factor(res.sig.plot$sporing)

# Plot all significant features (horizontal barchart)
print(ggplot()
      + geom_bar(data=res.sig.plot, aes(x=description_short, y=effect_size, fill=sporing), stat="identity", alpha=0.8, size=0.4)
      + geom_point(data=inf_df, aes(y = inf_vals, x = inf_idx), position = position_dodge(width = 0.9), shape=8, size=5, show.legend=FALSE)
      + theme_bw()
      + rremove("grid")
      + ggtitle(paste0(feat, "_", feat_type))
      + coord_flip()
      + theme(plot.title = element_text(hjust = 0.5))
      + ylab(y.lab)
      + scale_x_discrete(limits=res.sig.plot$description_short)
      + scale_fill_manual(values = c("red3", "dodgerblue3"))
      + theme(axis.title.y = element_blank()))

# Plot all significant features (vertical barchart)
res.sig.plot.rev <- res.sig.plot[order(res.sig.plot$effect_size, decreasing=FALSE),]
inf_df.rev <- inf_df
inf_df.rev$inf_idx <- nrow(res.sig.plot) - inf_df.rev$inf_idx + 1
inf_df.rev <- inf_df.rev[order(inf_df.rev$inf_idx, decreasing=FALSE),]

print(ggplot()
      + geom_bar(data=res.sig.plot.rev, aes(x=description_short, y=effect_size, fill=sporing), stat="identity", alpha=0.8, size=0.4)
      + geom_point(data=inf_df.rev, aes(y = inf_vals, x = inf_idx), position = position_dodge(width = 0.9), shape=8, size=5, show.legend=FALSE)
      + theme_bw()
      + rremove("grid")
      + ggtitle(paste0(feat, "_", feat_type))
      + ylab(y.lab)
      + scale_x_discrete(limits=res.sig.plot.rev$description_short)
      + scale_fill_manual(values = c("red3", "dodgerblue3"))
      + theme(plot.title = element_text(hjust = 0.5))
      + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))
