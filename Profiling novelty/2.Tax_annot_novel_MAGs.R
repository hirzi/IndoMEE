### This scripts adds taxonomic metadata to novel MAGs/new species identified from magscreen ###

# Load libraries
library(stringr)
library(gtools)

# Define assembly
assembly <- "hybrid"

# Read in genome metadata
if (assembly == "hybrid") {
  novel_MAGs <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/magscreen/Results/HybridAssembly_Comp50Cont5/new_species.txt")
  bin_tax_bac_raw <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/tax_annotation/taxonomy/GTDB_TK/HybridAssembly_Comp50Cont5/gtdbtk.bac120.summary.tsv", header = TRUE, sep = "\t")
  bin_tax_arc_raw <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/tax_annotation/taxonomy/GTDB_TK/HybridAssembly_Comp50Cont5/gtdbtk.ar53.summary.tsv", header = TRUE, sep = "\t")
  prefix <- "HybridAssembly_Comp50Cont5"
} else if (assembly == "co") {
  novel_MAGs <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/magscreen/Results/.../new_species.txt")
  bin_tax_bac_raw <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/tax_annotation/taxonomy/GTDB_TK/CoAssembly/gtdbtk.bac120.summary.tsv", header = TRUE, sep = "\t")
  bin_tax_arc_raw <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/tax_annotation/taxonomy/GTDB_TK/CoAssembly/gtdbtk.ar53.summary.tsv", header = TRUE, sep = "\t")
  prefix <- "CoAssembly_draft1"
} else if (assembly == "individual") {
  novel_MAGs <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/magscreen/Results/.../new_species.txt")
  bin_tax_bac_raw <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/tax_annotation/taxonomy/GTDB_TK/IndAssembly_draft2/gtdbtk.bac120.summary.tsv", header = TRUE, sep = "\t")
  bin_tax_arc_raw <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/tax_annotation/taxonomy/GTDB_TK/IndAssembly_draft2/gtdbtk.ar53.summary.tsv", header = TRUE, sep = "\t")
  prefix <- "IndAssembly_draft2"
} else if (assembly == "individual_sANI98") {
  novel_MAGs <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/magscreen/Results/.../new_species.txt")
  bin_tax_bac_raw <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/tax_annotation/taxonomy/GTDB_TK/IndAssembly_draft2_sANI98/gtdbtk.bac120.summary.tsv", header = TRUE, sep = "\t")
  bin_tax_arc_raw <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/tax_annotation/taxonomy/GTDB_TK/IndAssembly_draft2_sANI98/gtdbtk.ar53.summary.tsv", header = TRUE, sep = "\t")
  prefix <- "IndAssembly_draft2_sANI98"
}
bin_tax_raw <- rbind(bin_tax_bac_raw, bin_tax_arc_raw)

# Format dataframes
colnames(novel_MAGs) <- "genome"
novel_MAGs$genome=gsub(".fa","",novel_MAGs$genome)
colnames(bin_tax_raw)[1] <- "genome"
if (assembly == "co") {
  bin_tax_raw$genome=gsub("CO_","",bin_tax_raw$genome)
} else if (assembly == "individual" || assembly == "individual_sANI98") {
  bin_tax_raw$genome=gsub("IND_","",bin_tax_raw$genome)
}

# Merge dataframes
novel_MAGs_tax <- merge(novel_MAGs, bin_tax_raw, by = "genome")

# Write out metadata table and add genome metadata to metamap outputs
write.table(novel_MAGs_tax, paste0("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/magscreen/Results/", prefix, "/tax_new_species.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
