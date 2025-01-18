### This scripts adds metadata (inc. taxonomic IDs) to metamap outputs ###

# Load libraries
library(stringr)
library(gtools)

# Define assembly
assembly <- "hybrid"

# Read in genome metadata
if (assembly == "hybrid") {
  # bin_tax_order <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_ASSEMBLY/5.4.DREP/Comp50Cont5/dereplicated_genomes.txt", header = FALSE)
  # bin_stats_raw <- read.csv("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_ASSEMBLY/5.4.DREP/Comp50Cont5/data_tables/genomeInformation.csv", header = TRUE)
  # bin_tax_bac_raw <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/tax_annotation/taxonomy/GTDB_TK/HybridAssembly_Comp50Cont5/gtdbtk.bac120.summary.tsv", header = TRUE, sep = "\t")
  # bin_tax_arc_raw <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/tax_annotation/taxonomy/GTDB_TK/HybridAssembly_Comp50Cont5/gtdbtk.ar53.summary.tsv", header = TRUE, sep = "\t")
  # prefix <- "HybridAssembly_Comp50Cont5_draft1"
  bin_tax_order <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_ASSEMBLY/5.4.DREP/Comp50Cont5/dereplicated_genomes.txt", header = FALSE)
  bin_stats_raw <- read.csv("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_ASSEMBLY/5.4.DREP/Comp50Cont5/data_tables/genomeInformation.csv", header = TRUE)
  bin_tax_bac_raw <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/tax_annotation/taxonomy/GTDB_TK/HybridAssembly_sANI95_Comp50Cont5/gtdbtk.bac120.summary.tsv", header = TRUE, sep = "\t")
  bin_tax_arc_raw <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/tax_annotation/taxonomy/GTDB_TK/HybridAssembly_sANI95_Comp50Cont5/gtdbtk.ar53.summary.tsv", header = TRUE, sep = "\t")
  prefix <- "HybridAssembly_Comp50Cont5_draft2"
} else if (assembly == "hybrid_sANI98") {
  bin_tax_order <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_ASSEMBLY/5.4.DREP_sfastANI_98/Comp50Cont5/dereplicated_genomes.txt", header = FALSE)
  bin_stats_raw <- read.csv("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_ASSEMBLY/5.4.DREP_sfastANI_98/Comp50Cont5/data_tables/genomeInformation.csv", header = TRUE)
  bin_tax_bac_raw <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/tax_annotation/taxonomy/GTDB_TK/HybridAssembly_sANI98_Comp50Cont5/gtdbtk.bac120.summary.tsv", header = TRUE, sep = "\t")
  bin_tax_arc_raw <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/tax_annotation/taxonomy/GTDB_TK/HybridAssembly_sANI98_Comp50Cont5/gtdbtk.ar53.summary.tsv", header = TRUE, sep = "\t")
  prefix <- "HybridAssembly_sANI98_Comp50Cont5_draft2"
} else if (assembly == "co") {
  bin_tax_order <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/CO_ASSEMBLY/5.5.FILTER_GUNC/dereplicated_genomes_gunc_finalFiltered_names.txt", header = FALSE)
  bin_stats_raw <- read.csv("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/CO_ASSEMBLY/5.4.DREP/data_tables/genomeInformation.csv", header = TRUE)
  bin_tax_bac_raw <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/tax_annotation/taxonomy/GTDB_TK/CoAssembly/gtdbtk.bac120.summary.tsv", header = TRUE, sep = "\t")
  bin_tax_arc_raw <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/tax_annotation/taxonomy/GTDB_TK/CoAssembly/gtdbtk.ar53.summary.tsv", header = TRUE, sep = "\t")
  prefix <- "CoAssembly_draft1"
} else if (assembly == "co_sANI98") {
  bin_tax_order <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/CO_ASSEMBLY/5.5.FILTER_GUNC_sANI98/dereplicated_genomes_gunc_finalFiltered_names.txt", header = FALSE)
  bin_stats_raw <- read.csv("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/CO_ASSEMBLY/5.4.DREP_sANI98/data_tables/genomeInformation.csv", header = TRUE)
  bin_tax_bac_raw <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/tax_annotation/taxonomy/GTDB_TK/CoAssembly_sANI98/gtdbtk.bac120.summary.tsv", header = TRUE, sep = "\t")
  bin_tax_arc_raw <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/tax_annotation/taxonomy/GTDB_TK/CoAssembly_sANI98/gtdbtk.ar53.summary.tsv", header = TRUE, sep = "\t")
  prefix <- "CoAssembly_sANI98_draft2"
} else if (assembly == "individual") {
  bin_tax_order <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/5.5.FILTER_GUNC/dereplicated_genomes_gunc_finalFiltered_names.txt", header = FALSE) #renamed to 5.5.FILTER_GUNC_sANI95
  bin_stats_raw <- read.csv("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/5.4.DREP/data_tables/genomeInformation.csv", header = TRUE) #renamed to 5.4.DREP_sANI95
  bin_tax_bac_raw <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/tax_annotation/taxonomy/GTDB_TK/IndAssembly_draft2/gtdbtk.bac120.summary.tsv", header = TRUE, sep = "\t")
  bin_tax_arc_raw <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/tax_annotation/taxonomy/GTDB_TK/IndAssembly_draft2/gtdbtk.ar53.summary.tsv", header = TRUE, sep = "\t")
  prefix <- "IndAssembly_draft2"
} else if (assembly == "individual_sANI98") {
  bin_tax_order <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/5.5.FILTER_GUNC_sANI_98/dereplicated_genomes_gunc_finalFiltered_names.txt", header = FALSE)
  bin_stats_raw <- read.csv("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/5.4.DREP_sANI_98/data_tables/genomeInformation.csv", header = TRUE)
  bin_tax_bac_raw <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/tax_annotation/taxonomy/GTDB_TK/IndAssembly_draft2_sANI98/gtdbtk.bac120.summary.tsv", header = TRUE, sep = "\t")
  bin_tax_arc_raw <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/tax_annotation/taxonomy/GTDB_TK/IndAssembly_draft2_sANI98/gtdbtk.ar53.summary.tsv", header = TRUE, sep = "\t")
  prefix <- "IndAssembly_draft2_sANI98"
} else if (assembly == "local_test") {
  #bin_tax_order <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Script working files/dereplicated_genomes.txt", header = FALSE)
  #bin_stats_raw <- read.csv("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Script working files/genomeInformation.csv", header = TRUE)
  #bin_tax_bac_raw <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Script working files/gtdbtk.bac120.summary.tsv", header = TRUE, sep = "\t")
  #bin_tax_arc_raw <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Script working files/gtdbtk.ar53.summary.tsv", header = TRUE, sep = "\t")
  #prefix <- "HybridAssembly_Comp50Cont5"
  # bin_tax_order <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Script working files/CoAssembly_draft1/dereplicated_genomes_gunc_finalFiltered_names.txt", header = FALSE)
  # bin_stats_raw <- read.csv("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Script working files/CoAssembly_draft1/genomeInformation.csv", header = TRUE)
  # bin_tax_bac_raw <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Script working files/CoAssembly_draft1/gtdbtk.bac120.summary.tsv", header = TRUE, sep = "\t")
  # bin_tax_arc_raw <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Script working files/CoAssembly_draft1/gtdbtk.ar53.summary.tsv", header = TRUE, sep = "\t")
  # prefix <- "CoAssembly_draft1"
  # bin_tax_order <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Script working files/IndAssembly_draft2/dereplicated_genomes_gunc_finalFiltered_names.txt", header = FALSE)
  # bin_stats_raw <- read.csv("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Script working files/IndAssembly_draft2/genomeInformation.csv", header = TRUE)
  # bin_tax_bac_raw <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Script working files/IndAssembly_draft2/gtdbtk.bac120.summary.tsv", header = TRUE, sep = "\t")
  # bin_tax_arc_raw <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Script working files/IndAssembly_draft2/gtdbtk.ar53.summary.tsv", header = TRUE, sep = "\t")
  # prefix <- "IndAssembly_draft2"
}
bin_tax_raw <- rbind(bin_tax_bac_raw, bin_tax_arc_raw)

# Format dataframes
colnames(bin_tax_order) <- "genome"
bin_tax_order$genome=gsub(".fa","",bin_tax_order$genome)
colnames(bin_tax_raw)[1] <- "genome"
bin_stats_raw$genome=gsub(".fa","",bin_stats_raw$genome)
if (assembly == "co" || assembly == "co_sANI98") {
  bin_tax_raw$genome=gsub("CO_","",bin_tax_raw$genome)
} else if (assembly == "individual" || assembly == "individual_sANI98") {
  bin_tax_raw$genome=gsub("IND_","",bin_tax_raw$genome)
}

# Merge dataframes
genomes_selected_metadata_raw <- merge(bin_tax_raw, bin_stats_raw, by = "genome")

# Order dataframe according to input file used to make custom MAG reference database (i.e. in 99.Make_customMAGRef.sh, e.g. dereplicated_genomes.txt)
genomes_selected_metadata_raw <- genomes_selected_metadata_raw[match(bin_tax_order$genome, genomes_selected_metadata_raw$genome),]
rownames(genomes_selected_metadata_raw) <- seq(1,nrow(genomes_selected_metadata_raw))
# Check order is correct
all(genomes_selected_metadata_raw$genome == bin_tax_order$genome)

# Add final MAG identifiers (as in custom MAG reference database)
genomes_selected_metadata_raw$bin <- paste0("MGYG",sprintf("%09d", as.numeric(rownames(genomes_selected_metadata_raw))))
genomes_selected_metadata_raw <- genomes_selected_metadata_raw[,c(ncol(genomes_selected_metadata_raw),1:(ncol(genomes_selected_metadata_raw)-1))]
colnames(genomes_selected_metadata_raw)[c(1,2)] <- c("Genome", "original_bin")

# Write out metadata table and add genome metadata to metamap outputs
if (assembly == "local_test") {
  # write.table(genomes_selected_metadata_raw, paste0("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Script working files/", prefix,"/metawrap_", prefix, "_metadata.tsv"), sep = "\t", row.names = FALSE)
  # setwd(paste0("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Script working files/", prefix, "/summary"))
} else {
  write.table(genomes_selected_metadata_raw, paste0("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metamap/reference_assemblies_metadata/metawrap_", prefix, "_metadata.tsv"), sep = "\t", row.names = FALSE)
  setwd(paste0("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metamap/metamap_on_metaWrap_MAGs_results/", prefix, "/summary"))
}
metamap_output_files <- c("bwa_counts-total.csv", "bwa_counts-unique.csv", "bwa_cov-est.csv", "bwa_cov-exp.csv")
#f <- "bwa_counts-total.csv"
for (f in metamap_output_files) {
  # Read in metamap output
  out_raw <- read.csv(f)
  # Remove .fa from genome names
  out_raw$Genome <- gsub(".fa", "", out_raw$Genome)
  # Merge metamap outputs with metadata
  out <- merge(genomes_selected_metadata_raw, out_raw, by = "Genome")
  # Sort by total abundance
  out$SumAcrossAllSamples <- rowSums(out[,(ncol(genomes_selected_metadata_raw)+1):ncol(out)])
  out <- out[order(out$SumAcrossAllSamples, decreasing = TRUE),]
  # Write out tables
  write.table(out, paste0(str_split(f, ".csv")[[1]][1], "_wMetadata.tsv"), sep = "\t", row.names = FALSE)
}
