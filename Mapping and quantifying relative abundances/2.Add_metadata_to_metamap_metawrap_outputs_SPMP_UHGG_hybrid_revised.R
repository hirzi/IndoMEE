### This scripts adds metadata (inc. taxonomic IDs) to metamap outputs ###

# Load libraries
library(stringr)
library(gtools)

# Define assembly
assembly <- "SPMP_sANI95"

# Read in genome metadata
if (assembly == "SPMP" || assembly == "SPMP_sANI95" || assembly == "hybrid_UHGG_SPMP") {
  #genomes_selected_metadata_raw_1 <-read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assemblies/SPMP/MAGs_tax_metadata_41467_2022_33782_MOESM9_ESM.txt", header = TRUE, sep = "\t")
  genomes_selected_metadata_raw_1 <-read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Public_genDiversity_databases/SPMP/MAGs_tax_metadata_41467_2022_33782_MOESM9_ESM.txt", header = TRUE, sep = "\t")
  #genomes_selected_metadata_raw_2 <-read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assemblies/SPMP/SPMPPolishedMAGs.txt", header = FALSE, sep = "\t")
  genomes_selected_metadata_raw_2 <-read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Public_genDiversity_databases/SPMP/SPMPPolishedMAGs.txt", header = FALSE, sep = "\t")
  genomes_selected_metadata_raw_2[,2] <- paste0("MGYG",sprintf("%09d", as.numeric(rownames(genomes_selected_metadata_raw_2))))
  #genomes_selected_metadata_raw_3 <-read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assemblies/SPMP/SPMPPolishedMAGs_lengths.txt", header = FALSE, sep = "\t")
  genomes_selected_metadata_raw_3 <-read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Public_genDiversity_databases/SPMP/SPMPPolishedMAGs_lengths.txt", header = FALSE, sep = "\t")
  colnames(genomes_selected_metadata_raw_3) <- c("MAG.ID", "Length")
  colnames(genomes_selected_metadata_raw_2) <- c("MAG.ID", "Genome")
  genomes_selected_metadata_raw_2$MAG.ID <- gsub(".fasta", "", genomes_selected_metadata_raw_2$MAG.ID)
  genomes_selected_metadata_raw_3$MAG.ID <- gsub(".fasta", "", genomes_selected_metadata_raw_3$MAG.ID)
  genomes_selected_metadata_raw_1$MAG.ID <- gsub("SPMP", "TLL", genomes_selected_metadata_raw_1$MAG.ID)
  genomes_selected_metadata_raw <- merge(genomes_selected_metadata_raw_1, genomes_selected_metadata_raw_2, by = "MAG.ID")
  genomes_selected_metadata_raw <- merge(genomes_selected_metadata_raw, genomes_selected_metadata_raw_3, by = "MAG.ID")
  genomes_selected <- genomes_selected_metadata_raw[,c(ncol(genomes_selected_metadata_raw) - 1 , 1:(ncol(genomes_selected_metadata_raw)-2), ncol(genomes_selected_metadata_raw))]
}

if (assembly == "SPMP" || assembly == "SPMP_sANI95") {
  prefix <- assembly
} else if (assembly == "hybrid_UHGG") {
  #metadata_UHGG <-read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/UHGG database/Gut/genomes-selected_metadata.tsv", header = TRUE)
  metadata_UHGG <-read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Microbial_references/MGnifyGenomes_humanGut/v2.0.1/genomes-selected_metadata.tsv", header = TRUE)
  #metadata_hybrid <-read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assemblies/HybridAssembly_Comp50Cont5_draft1/metawrap_HybridAssembly_Comp50Cont5_metadata.tsv", header = TRUE)
  #metadata_hybrid <-read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metamap/reference_assemblies_metadata/metawrap_HybridAssembly_Comp50Cont5_draft1_metadata.tsv", header = TRUE)
  metadata_hybrid <-read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metamap/reference_assemblies_metadata/metawrap_HybridAssembly_Comp50Cont5_draft2_metadata.tsv", header = TRUE)
  #genomes <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assemblies/Hybrid_UHGG_Assembly_Comp50Cont5_draft1/dereplicated_genomes.txt", header = FALSE)
  genomes <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_UHGG_ASSEMBLY/5.4.DREP/dereplicated_genomes.txt", header = FALSE)
  prefix <- "Hybrid_UHGG_Assembly_Comp50Cont5_draft2"
} else if (assembly == "hybrid_UHGG_SPMP") {
  #metadata_UHGG <-read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/UHGG database/Gut/genomes-selected_metadata.tsv", header = TRUE)
  metadata_UHGG <-read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Microbial_references/MGnifyGenomes_humanGut/v2.0.1/genomes-selected_metadata.tsv", header = TRUE)
  #metadata_hybrid <-read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assemblies/HybridAssembly_Comp50Cont5_draft1/metawrap_HybridAssembly_Comp50Cont5_metadata.tsv", header = TRUE)
  #metadata_hybrid <-read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metamap/reference_assemblies_metadata/metawrap_HybridAssembly_Comp50Cont5_draft1_metadata.tsv", header = TRUE)
  metadata_hybrid <-read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metamap/reference_assemblies_metadata/metawrap_HybridAssembly_Comp50Cont5_draft2_metadata.tsv", header = TRUE)
  metadata_SPMP <- genomes_selected
  #genomes <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assemblies/Hybrid_UHGG_SPMP_Assembly_Comp50Cont5_draft1/dereplicated_genomes.txt", header = FALSE)
  genomes <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_UHGG_SPMP_ASSEMBLY/5.4.DREP/dereplicated_genomes.txt", header = FALSE)
  prefix <- "Hybrid_UHGG_SPMP_Assembly_Comp50Cont5_draft2"
}

if (assembly == "hybrid_UHGG" || assembly == "hybrid_UHGG_SPMP") {
  genomes[,2] <- paste0("MGYG",sprintf("%09d", as.numeric(rownames(genomes))))
  colnames(genomes) <- c("Genome", "Final_MAG_ID")
  genomes$Genome <- gsub(".fasta", ".fa", genomes$Genome)
  genomes$Genome <- gsub(".fa", "", genomes$Genome)
  colnames(metadata_hybrid)[c(1,2)] <- c("Original_MAG_ID", "Genome")
  genomes_selected_hybrid <- merge(genomes, metadata_hybrid, by = "Genome")
  genomes_selected_UHGG <- merge(genomes, metadata_UHGG, by = "Genome")
  UHGG_selCols <- c("Genome", "Final_MAG_ID", "Species_rep", "Lineage", "Genome_type", "Length", "N50", "Completeness", "Contamination")
  genomes_selected_UHGG <- genomes_selected_UHGG[ , (names(genomes_selected_UHGG) %in% UHGG_selCols)]
  genomes_selected_UHGG <- genomes_selected_UHGG[,c(2, 8, 1, 9, 3, 4, 5, 6, 7)]
  colnames(genomes_selected_UHGG) <- c("Final_MAG_ID", "Original_MAG_ID", "Original_bin_ID", "Taxonomy", "Source", "Length", "N50", "Completeness", "Contamination")
  hybrid_selCols <- c("Genome", "Final_MAG_ID", "Original_MAG_ID", "classification", "Genome_type", "length", "N50", "completeness", "contamination")
  genomes_selected_hybrid <- genomes_selected_hybrid[ , (names(genomes_selected_hybrid) %in% hybrid_selCols)]
  genomes_selected_hybrid[,ncol(genomes_selected_hybrid) + 1] <- "MAG"
  genomes_selected_hybrid <- genomes_selected_hybrid[,c(2, 3, 1, 4, 9, 7, 8, 5, 6)]
  colnames(genomes_selected_hybrid) <- c("Final_MAG_ID", "Original_MAG_ID", "Original_bin_ID", "Taxonomy", "Source", "Length", "N50", "Completeness", "Contamination")
  genomes_selected <- rbind(genomes_selected_hybrid, genomes_selected_UHGG)
}

if (assembly == "hybrid_UHGG_SPMP") {
  genomes_selected_hybrid_UHGG <- genomes_selected
  colnames(metadata_SPMP)[c(1,2)] <- c("Original_MAG_ID", "Genome")
  genomes_selected_SPMP <- merge(genomes, metadata_SPMP, by = "Genome")
  genomes_selected_SPMP[, ncol(genomes_selected_SPMP) + 1 ] <- paste0(genomes_selected_SPMP$GTDB.Genus, "; ", genomes_selected_SPMP$GTDB.Species)
  colnames(genomes_selected_SPMP)[ncol(genomes_selected_SPMP)] <- "Taxonomy"
  SPMP_selCols <- c("Genome", "Final_MAG_ID", "Original_MAG_ID", "Length", "Taxonomy")
  genomes_selected_SPMP <- genomes_selected_SPMP[ , (names(genomes_selected_SPMP) %in% SPMP_selCols)]
  colnames(genomes_selected_SPMP) <- c("Original_bin_ID", "Final_MAG_ID", "Original_MAG_ID", "Length", "Taxonomy")
  genomes_selected_SPMP <- genomes_selected_SPMP[,c(2, 3, 1, 5, 4)]
  hybrid_UHGG_selCols <- c("Original_bin_ID", "Final_MAG_ID", "Original_MAG_ID", "Length", "Taxonomy")
  genomes_selected_hybrid_UHGG <- genomes_selected_hybrid_UHGG[ , (names(genomes_selected_hybrid_UHGG) %in% hybrid_UHGG_selCols)]
  genomes_selected <- rbind(genomes_selected_hybrid_UHGG, genomes_selected_SPMP)
}

if (assembly == "hybrid_UHGG" || assembly == "hybrid_UHGG_SPMP") {
  colnames(genomes_selected)[c(1)] <- "Genome"
}

# Write out metadata table and add genome metadata to metamap outputs
#write.table(genomes_selected, paste0("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assemblies/", prefix,"/metawrap_", prefix, "_metadata.tsv"), sep = "\t", row.names = FALSE)
#setwd(paste0("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assemblies/", prefix, "/summary"))
write.table(genomes_selected, paste0("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metamap/reference_assemblies_metadata/metawrap_", prefix, "_metadata.tsv"), sep = "\t", row.names = FALSE)
setwd(paste0("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metamap/metamap_on_metaWrap_MAGs_results/", prefix, "/summary"))

metamap_output_files <- c("bwa_counts-total.csv", "bwa_counts-unique.csv", "bwa_cov-est.csv", "bwa_cov-exp.csv")
#f <- "bwa_counts-total.csv"
for (f in metamap_output_files) {
  # Read in metamap output
  out_raw <- read.csv(f)
  # Remove .fa from genome names
  out_raw$Genome <- gsub(".fa", "", out_raw$Genome)
  # Merge metamap outputs with metadata
  out <- merge(genomes_selected, out_raw, by = "Genome")
  # Sort by total abundance
  out$SumAcrossAllSamples <- rowSums(out[,(ncol(genomes_selected)+1):ncol(out)])
  out <- out[order(out$SumAcrossAllSamples, decreasing = TRUE),]
  # Write out tables
  write.table(out, paste0(str_split(f, ".csv")[[1]][1], "_wMetadata.tsv"), sep = "\t", row.names = FALSE)
}
