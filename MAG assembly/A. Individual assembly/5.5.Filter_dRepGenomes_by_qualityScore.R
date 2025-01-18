### This script filters a list of genomes based on their completeness and contamination rates ###

# Read in list of dereplicated genomes and checkm quality scores
checkm_df_raw <- read.csv("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/5.2.CHECKM2/checkm2_genInfo.csv")
derep_genomes <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/5.5.FILTER_GUNC/dereplicated_genomes.txt")
colnames(derep_genomes) <- "genome"

# Filter genomes on S = Completion - 5*Contamination >= 50%
checkm_df_raw_filtered <- merge(derep_genomes, checkm_df_raw, by = "genome")
checkm_df_raw_filtered[,ncol(checkm_df_raw_filtered) + 1] <- checkm_df_raw_filtered$completeness - (5 * checkm_df_raw_filtered$contamination)
colnames(checkm_df_raw_filtered)[ncol(checkm_df_raw_filtered)] <- "quality_score"
checkm_df_raw_filtered2 <- checkm_df_raw_filtered[checkm_df_raw_filtered$quality_score >= 50,]

# Write out filtered table
write.table(checkm_df_raw_filtered2, "/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/5.5.FILTER_GUNC/dereplicated_genomes_filtered.txt", row.names = FALSE, quote = FALSE, sep = "\t")
