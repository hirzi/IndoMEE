### This script filters out gunc-identified chimeric genomes from a list of genomes  ###

# Load libraries
library(stringr)

# Read in drep cluster tabnle, checkm scores and gunc scores 
checkm_raw <- read.csv("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/5.2.CHECKM2/checkm2_genInfo.csv")
drep_clusters_raw <- read.csv("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/5.4.DREP/data_tables/Cdb.csv")
derep_genomes_filtered <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/5.5.FILTER_GUNC/dereplicated_genomes_filtered.txt", header = TRUE)
gunc_raw <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/5.5.FILTER_GUNC/GUNC.progenomes_2.1.maxCSS_level.tsv", header = TRUE)

# Get number of genomes represented in each secondary cluster.
# First, split secondary cluster column 
drep_clusters_raw_split <- cbind(drep_clusters_raw, str_split_fixed(drep_clusters_raw$secondary_cluster, "_", 2))
colnames(drep_clusters_raw_split)[c(ncol(drep_clusters_raw_split)-1, ncol(drep_clusters_raw_split))] <- c("sec_cluster_prefix", "sec_cluster_suffix")
drep_clusters_raw_split$sec_cluster_prefix <- as.numeric(drep_clusters_raw_split$sec_cluster_prefix)
# Then enumerate number of genomes for each secondary cluster. For reference, see: https://stackoverflow.com/questions/57125190/remove-singleton-entry-from-table-in-r 
TAB=table(drep_clusters_raw_split$sec_cluster_prefix)
multitons <- drep_clusters_raw_split[ifelse(TAB[drep_clusters_raw_split$sec_cluster_prefix]==1, FALSE, TRUE),]
multitons[,ncol(multitons) + 1] <- "multiton"
colnames(multitons)[ncol(multitons)] <- "frequency"  
singletons <- drep_clusters_raw_split[ifelse(TAB[drep_clusters_raw_split$sec_cluster_prefix]==1, TRUE, FALSE),]
singletons[,ncol(singletons) + 1] <- "singleton"
colnames(singletons)[ncol(singletons)] <- "frequency" 
drep_clusters_raw_splitFreq <- rbind(multitons, singletons)

# Merge dataframes
df_merged_raw <- merge(derep_genomes_filtered, drep_clusters_raw_splitFreq, by = "genome")
df_merged_raw <- cbind(df_merged_raw, str_split_fixed(df_merged_raw$genome, ".fa", 2)[,1])
colnames(df_merged_raw)[1] <- "name"
colnames(df_merged_raw)[ncol(df_merged_raw)] <- "genome"
df_merged_raw <- merge(df_merged_raw, gunc_raw, by = "genome")

# Filter out genomes that match all of the following criteria: i) gunc flagged, ii) are singletons (i.e. dRep clusters with only one member) and ii) < 50% or 90% completion. For reference, see: https://github.com/alexmsalmeida/magscreen)
df_merged_filtered <- df_merged_raw[df_merged_raw$pass.GUNC == "False",]
df_merged_filtered <- df_merged_filtered[df_merged_filtered$frequency == "singleton",]
df_merged_filtered <- df_merged_filtered[df_merged_filtered$completeness < 90,]

# Get the complement of those filtered out to get a list of genomes we retain
df_merged_kept <- df_merged_raw[!df_merged_raw$genome %in% df_merged_filtered$genome, ]

# Write out filtered table
write.table(df_merged_kept, "/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/5.5.FILTER_GUNC/dereplicated_genomes_gunc_finalFiltered.txt", row.names = FALSE, quote = FALSE, sep = "\t")

