# Load libraries
library(data.table)
library(matrixStats)
library(dplyr)
library(purrr)
library(ggplot2)

# Load data - proteins
setwd("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/Functional enrichment/")
eggnog.raw = read.delim("HybridAssembly_Comp50Cont5_Summary/eggnog_emapper_summary.tsv", header=FALSE, stringsAsFactors = FALSE)[,c(1,7)][-1,]
colnames(eggnog.raw) <- c("seq_id", "cog")
eggnog.proteins <- eggnog.raw$seq_id
#cog.descr = read.delim("prokaryotes/functions/cog_descriptions.tsv", stringsAsFactors = FALSE)
all.proteins <- scan("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/Functional enrichment/HybridAssembly95_vs_UHGGallMAGs_allMAGs_allProteins", what = "")
all.proteins.prefix <-  sub("_[^_]*$", "", all.proteins)
novel.proteins <- all.proteins[!(all.proteins %in% eggnog.proteins)]
novel.proteins.prefix <-  sub("_[^_]*$", "", novel.proteins)
unchar.proteins <- eggnog.raw[eggnog.raw$cog == "S", ]$seq_id
unchar.proteins.prefix <-  sub("_[^_]*$", "", unchar.proteins)

# Load data - MAGs
magscreen_dir <- "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/magscreen/Novel_MAGs/"
magscreen_comp <- "HybridAssembly_sANI95_Comp50Cont5_vs_UHGG_allMAGs"
#magscreen_comp <- "HybridAssembly_sANI98_Comp50Cont5_vs_UHGG_allMAGs"
novel.MAGs <- read.table(paste0(magscreen_dir, magscreen_comp, "/new_species.txt"))$V1
novel.MAGs <- gsub(".fa", "", novel.MAGs)

# Partition into novel and known MAG sets
n.novel.MAGs <- length(novel.MAGs)
n.all.MAGs <- 1304
n.known.MAGs <- n.all.MAGs - n.novel.MAGs
all.proteins.novel.MAGs <- all.proteins[(all.proteins.prefix %in% novel.MAGs)]
all.proteins.known.MAGs <- all.proteins[!(all.proteins.prefix %in% novel.MAGs)]
novel.proteins.novel.MAGs <- novel.proteins[(novel.proteins.prefix %in% novel.MAGs)]
novel.proteins.known.MAGs <- novel.proteins[!(novel.proteins.prefix %in% novel.MAGs)]
unchar.proteins.novel.MAGs <- unchar.proteins[(unchar.proteins.prefix %in% novel.MAGs)]
unchar.proteins.known.MAGs <- unchar.proteins[!(unchar.proteins.prefix %in% novel.MAGs)]


# Normalise
length(novel.proteins.novel.MAGs) / length(all.proteins.novel.MAGs) # % proteins that are novel in novel MAG set
length(novel.proteins.known.MAGs) / length(all.proteins.known.MAGs) # % proteins that are novel in known MAG set
length(unchar.proteins.novel.MAGs) / length(all.proteins.novel.MAGs) # % proteins that are uncharacterised in novel MAG set
length(unchar.proteins.known.MAGs) / length(all.proteins.known.MAGs) # % proteins that are uncharacterised in known MAG set

# Fisher's exact or chi-square test
cont_table_novelProteins <- matrix(c(length(novel.proteins.novel.MAGs), length(novel.proteins.known.MAGs), (length(all.proteins.novel.MAGs) - length(novel.proteins.novel.MAGs)), (length(all.proteins.known.MAGs) - length(novel.proteins.known.MAGs))), nrow = 2,
                     dimnames = list(c("Novel MAG", "Known MAG"), c("Novel proteins", "Known proteins")))
cont_table_uncharProteins <- matrix(c(length(unchar.proteins.novel.MAGs), length(unchar.proteins.known.MAGs), (length(all.proteins.novel.MAGs) - length(unchar.proteins.novel.MAGs) - length(novel.proteins.novel.MAGs)), (length(all.proteins.known.MAGs) - length(unchar.proteins.known.MAGs) - length(novel.proteins.known.MAGs))), nrow = 2,
                                   dimnames = list(c("Novel MAG", "Known MAG"), c("Uncharacterised proteins", "Characterised proteins")))
fisher.test(cont_table_novelProteins)
fisher.test(cont_table_uncharProteins)
chisq.test(cont_table_novelProteins)
chisq.test(cont_table_uncharProteins)

