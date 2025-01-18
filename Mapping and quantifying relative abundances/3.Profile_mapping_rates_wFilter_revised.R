### This script calculates ...
## For reference, see: https://github.com/alexmsalmeida/metamap ; https://instrain.readthedocs.io/en/latest/important_concepts.html#detecting-organisms-in-metagenomic-data

# Define whether metamap was run on the UHGG database or on your metawrap-assembled MAGs
cpu <- "remote"
#assembly_list <- c("UHGG", "IndAssembly_draft2", "IndAssembly_draft2_sANI98", "CoAssembly_draft1", "CoAssembly_sANI98_draft2", "HybridAssembly_Comp50Cont5_draft2", "HybridAssembly_sANI98_Comp50Cont5_draft2", "SPMP", "Hybrid_UHGG_Assembly_Comp50Cont5_draft2", "Hybrid_UHGG_SPMP_Assembly_Comp50Cont5_draft2")
#assembly_list <- c("CoAssembly_sANI98_draft2", "HybridAssembly_Comp50Cont5_draft2", "HybridAssembly_sANI98_Comp50Cont5_draft2", "Hybrid_UHGG_Assembly_Comp50Cont5_draft2", "Hybrid_UHGG_SPMP_Assembly_Comp50Cont5_draft2")
assembly_list <- c("UHGG", "SPMP", "HybridAssembly_Comp50Cont5_draft2", "Hybrid_UHGG_SPMP_Assembly_Comp50Cont5_draft2")
#assembly <- "Hybrid_UHGG_Assembly_Comp50Cont5_draft2"
filter <- "weak"

# Define directory prefixes
if (cpu == "local") {
  results_dir <- "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assemblies/On_hostDecontaminated_fastqs/"
  metadata_dir <- "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Metadata/"
} else if (cpu == "remote") {
  results_dir <- "/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metamap/metamap_on_metaWrap_MAGs_results/"
  metadata_dir <- "/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metamap/reference_assemblies_metadata/"
}

# For loop across assemblies
for (assembly in assembly_list) {
  
  prefix <- assembly
  print(prefix)
  
  # Set working directory
  setwd(paste0(results_dir, prefix, "/summary/"))
  #setwd(paste0(results_dir, prefix, "/summary_wSPMP/"))
  genomes_selected_metadata_raw <-read.table(paste0(metadata_dir, prefix, "/metawrap_", prefix, "_metadata.tsv"), header = TRUE)
  
  # Read in metamap outputs
  bwa_counts_total <- read.csv("bwa_counts-total.csv")   # read counts per genome across all samples including multi-mapped reads.
  bwa_counts_unique <- read.csv("bwa_counts-unique.csv")   # read counts per genome across all samples excluding multi-mapped reads.
  bwa_cov_est <- read.csv("bwa_cov-est.csv")   # breadth of coverage per genome across all samples.
  bwa_cov_exp <- read.csv("bwa_cov-exp.csv")   # expected breadth of coverage per genome across all samples based on their level of read depth.
  # Remove .fa from genome names
  if(assembly == "UHGG") {
    bwa_counts_total$Genome <- gsub(".fa", "", bwa_counts_total$Genome)
    bwa_counts_unique$Genome <- gsub(".fa", "", bwa_counts_unique$Genome)
    bwa_cov_est$Genome <- gsub(".fa", "", bwa_cov_est$Genome)
    bwa_cov_exp$Genome <- gsub(".fa", "", bwa_cov_exp$Genome)
  }
  # Remove X from sample names
  colnames(bwa_counts_total) <- gsub("X", "", colnames(bwa_counts_total))
  colnames(bwa_counts_unique) <- gsub("X", "", colnames(bwa_counts_unique))
  colnames(bwa_cov_est) <- gsub("X", "", colnames(bwa_cov_est))
  colnames(bwa_cov_exp) <- gsub("X", "", colnames(bwa_cov_exp))
  # Sort columns by sample name
  bwa_counts_total <- bwa_counts_total[ , c(1, (order(colnames(bwa_counts_total)[c(2:ncol(bwa_counts_total))]) + 1))]
  bwa_counts_unique <- bwa_counts_unique[ , c(1, (order(colnames(bwa_counts_unique)[c(2:ncol(bwa_counts_unique))]) + 1))]
  bwa_cov_est <- bwa_cov_est[ , c(1, (order(colnames(bwa_cov_est)[c(2:ncol(bwa_cov_est))]) + 1))]
  bwa_cov_exp <- bwa_cov_exp[ , c(1, (order(colnames(bwa_cov_exp)[c(2:ncol(bwa_cov_exp))]) + 1))]
  
  # Initialise filtered dfs
  bwa_counts_unique_filtered <- bwa_counts_unique
  bwa_counts_total_filtered <- bwa_counts_total
  
  # Check that columns and rows are ordered identically
  stopifnot(all(colnames(bwa_cov_est) == colnames(bwa_cov_exp)) == all(colnames(bwa_counts_unique) == colnames(bwa_counts_total)))
  stopifnot(all(rownames(bwa_cov_est) == rownames(bwa_cov_exp)) == all(rownames(bwa_counts_unique) == rownames(bwa_counts_total)))
  
  # Because of mis-mapping, abundance or counts > 0 doesn't necessarily imply the species is present (see: https://instrain.readthedocs.io/en/latest/important_concepts.html#detecting-organisms-in-metagenomic-data).
  # So we filter metamap counts outputs. Here we use a 50% breadth threshold (or 25% breadth threshold & >= 0.8 ratio estimated breadth: expected breadth) to determine presence/absence instead.
  # Alternatively, Alex suggested using the following threshold: > 5% breadth and 0.3 ratio of estimated breadth: expected breadth
  for (i in (seq(1, nrow(bwa_cov_est)))) {
    for (j in seq(2, ncol(bwa_cov_est))) {
      if(filter == "strong") {
        if((bwa_cov_est[i,j] >= 50) || (((bwa_cov_est[i,j] / bwa_cov_exp[i,j]) >= 0.8) && (bwa_cov_est[i,j] >= 25))) {
          #bwa_counts_unique_filtered[i,j] <- bwa_counts_unique_filtered[i,j]
          #bwa_counts_total_filtered[i,j] <- bwa_counts_total_filtered[i,j]
        } else {
          bwa_counts_unique_filtered[i,j] <- 0
          bwa_counts_total_filtered[i,j] <- 0
        }
      } else if (filter == "weak") {
        if((((bwa_cov_est[i,j] / bwa_cov_exp[i,j]) >= 0.3) && (bwa_cov_est[i,j] >= 5))) {
          #bwa_counts_unique_filtered[i,j] <- bwa_counts_unique_filtered[i,j]
          #bwa_counts_total_filtered[i,j] <- bwa_counts_total_filtered[i,j]
        } else {
          bwa_counts_unique_filtered[i,j] <- 0
          bwa_counts_total_filtered[i,j] <- 0
        }
      }
    }
  }
  
  if (assembly != "NAN") {
    # Merge filtered metamap outputs with metadata
    bwa_counts_unique_filtered_wMetadata <- merge(genomes_selected_metadata_raw, bwa_counts_unique_filtered, by = "Genome")
    bwa_counts_total_filtered_wMetadata <- merge(genomes_selected_metadata_raw, bwa_counts_total_filtered, by = "Genome")
    # Sort by total abundance
    bwa_counts_unique_filtered_wMetadata$SumAcrossAllSamples <- rowSums(bwa_counts_unique_filtered_wMetadata[,(ncol(genomes_selected_metadata_raw)+1):ncol(bwa_counts_unique_filtered_wMetadata)])
    bwa_counts_total_filtered_wMetadata$SumAcrossAllSamples <- rowSums(bwa_counts_total_filtered_wMetadata[,(ncol(genomes_selected_metadata_raw)+1):ncol(bwa_counts_total_filtered_wMetadata)])
  } else if (assembly == "NAN") { # i.e. if no metadata is available for assembly
    # Merge filtered metamap outputs with metadata
    bwa_counts_unique_filtered_wMetadata <- bwa_counts_unique_filtered
    bwa_counts_total_filtered_wMetadata <- bwa_counts_total_filtered
    # Sort by total abundance
    bwa_counts_unique_filtered_wMetadata$SumAcrossAllSamples <- rowSums(bwa_counts_unique_filtered_wMetadata[,2:ncol(bwa_counts_unique_filtered_wMetadata)])
    bwa_counts_total_filtered_wMetadata$SumAcrossAllSamples <- rowSums(bwa_counts_total_filtered_wMetadata[,2:ncol(bwa_counts_total_filtered_wMetadata)])
  }
  bwa_counts_unique_filtered_wMetadata <- bwa_counts_unique_filtered_wMetadata[order(bwa_counts_unique_filtered_wMetadata$SumAcrossAllSamples, decreasing = TRUE),]
  bwa_counts_total_filtered_wMetadata <- bwa_counts_total_filtered_wMetadata[order(bwa_counts_total_filtered_wMetadata$SumAcrossAllSamples, decreasing = TRUE),]
  
  # Count number of species
  species_richness_unique <- as.data.frame(colSums(bwa_counts_unique_filtered[,c(2:ncol(bwa_counts_unique_filtered))] != 0))
  species_richness_total <- as.data.frame(colSums(bwa_counts_total_filtered[,c(2:ncol(bwa_counts_total_filtered))] != 0))
  colnames(species_richness_unique)[ncol(species_richness_unique)] <- "nSpecies"
  colnames(species_richness_total)[ncol(species_richness_total)] <- "nSpecies"
  species_richness_unique$sample <- rownames(species_richness_unique)
  species_richness_total$sample <- rownames(species_richness_total)
  
  # Write out tables
  write.table(bwa_counts_unique_filtered_wMetadata, paste0("bwa_counts_unique_filtered_wMetadata_", prefix, "_", filter, "Filter.tsv"), sep = "\t", row.names = FALSE)
  write.table(bwa_counts_total_filtered_wMetadata, paste0("bwa_counts_total_filtered_wMetadata_", prefix, "_", filter, "Filter.tsv"), sep = "\t", row.names = FALSE)
  
  #### Let's compare with metaphlan (i.e. in terms of estimated species richness)
  # Initialise empty lists to hold results
  species_richness_sacParameter_list <- list()
  sac_fitdist_params_distmat_list <- list()
  sac_fitdist_params_richness_distmat_list <- list()
  
  # Import SGB (or other taxon level) abundance data
  if (cpu == "local") {
    microbiota_composition_data_raw <- read.csv("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Biobakery/Metaphlan_combined_pilotSamples_wTaxoInfo_rawTable.csv")
  } else if (cpu == "remote") {
    microbiota_composition_data_raw <- read.csv("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metamap/Metaphlan_combined_pilotSamples_wTaxoInfo_rawTable.csv")
  }
  nTaxonLevel <- 8
  nSGB <- ncol(microbiota_composition_data_raw) - nTaxonLevel
  # Filter samples
  samples_to_remove <- c("BAPDWEBECTRL", "CUST_POSCTR01", "CUST_POSCTR02", "NTCEBBATCH1", "NTCEB_BATCH2", "NTCEBBATCH3", "NTCNFW1", "PBCHFC001", "POSMOCKCTRL01", "POSMOCKCTRL02", "MALHULDG01", "MALHULDG03", "MALHULDG04", "PBSDGF001", "PBSDGF003", "PBSDGF005B", "PBSDGF005", "X124", "X130", "X131", "X238", "X239", "X367", "X379", "X388", "X391", "X394", "JKTCTRF01B2B3", "JKTCTRF01", "JKTCTRF04", "JKTCTRM02", "JKTCTRM03", "PBSHF001", "PBSHF002", "PBSHF003")
  microbiota_composition_data <- microbiota_composition_data_raw[, !(names(microbiota_composition_data_raw) %in% samples_to_remove)]
  
  # Filter taxa by abundance
  # For proper estimation of full SAC, you probably don't want to artificially truncate i.e. filter on abundance
  microbiota_composition_data <- microbiota_composition_data[,-c(2:8)]
  abundance_filter <- FALSE
  if (abundance_filter == TRUE) {
    abundance_filter_threshold <- 0.1
    microbiota_composition_data[, ncol(microbiota_composition_data)+1] <- rowSums(microbiota_composition_data[-c(1)])
    microbiota_composition_data[, ncol(microbiota_composition_data)+1] <- microbiota_composition_data[, ncol(microbiota_composition_data)]/(nSGB)
    colnames(microbiota_composition_data)[(ncol(microbiota_composition_data)-1):ncol(microbiota_composition_data)] <- c("Abundance", "NormAbundance")
    microbiota_composition_data_filtered <- microbiota_composition_data[microbiota_composition_data$NormAbundance >= abundance_filter_threshold, ]
    microbiota_composition_data_filtered <- microbiota_composition_data_filtered[,-c((ncol(microbiota_composition_data_filtered)-1):ncol(microbiota_composition_data_filtered))]
  } else if (abundance_filter == FALSE) {
    microbiota_composition_data_filtered <- microbiota_composition_data
  }
  
  # Calculate species richness
  species_richness_vec <- c()
  i <- 2
  for (ind in colnames(microbiota_composition_data_filtered)[2:ncol(microbiota_composition_data_filtered)]) {
    sac_species <- microbiota_composition_data_filtered[,c(1,i)]
    sac_species <- sac_species[sac_species[,2] != 0, ]
    species_richness_vec <- c(species_richness_vec, nrow(sac_species))
    i <- i + 1
  }
  species_richness_df <- as.data.frame(cbind(colnames(microbiota_composition_data_filtered)[2:ncol(microbiota_composition_data_filtered)], species_richness_vec))
  species_richness_df$species_richness_vec <- as.numeric(species_richness_df$species_richness_vec)
  colnames(species_richness_df) <- c("sample", "nSpecies")
  
  # Merge metaphlan and metamap estimates of species richness
  species_richness_summary <- merge(species_richness_df, species_richness_total, by = "sample")
  species_richness_summary <- merge(species_richness_summary, species_richness_unique, by = "sample")
  colnames(species_richness_summary) <- c("sample", "metaphlan", "metamap_total", "metamap_unique")
  
  # # And plot
  # plot(species_richness_summary$metaphlan, species_richness_summary$metamap_unique, pch=19)
  # #plot(species_richness_summary$metaphlan, species_richness_summary$metamap_total)
  
  # Let's calculate the percentage of unassigned reads
  # Note that here, unassigned reads would include human (host) reads, which we have not a priori mapped against (it might be useful to first map and remove human-mapped reads before running metamap, to more accurately quantify % dark matter. Alternatively we can get this figure from our metaphlan results).
  if (cpu == "local") {
    raw_reads <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/sample_nReads.txt")
    cleaned_reads <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/sample_nReads_cleanedFastQs.list")
  } else if (cpu == "remote") {
    raw_reads <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metamap/sample_nReads.txt")
    cleaned_reads <- read.table("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metamap/sample_nReads_cleanedFastQs.list")
  }
  colnames(raw_reads) <- c("sample", "fastq_raw_reads")
  colnames(cleaned_reads) <- c("sample", "fastq_cleaned_reads")
  raw_reads$fastq_raw_reads <- as.numeric(as.character(raw_reads$fastq_raw_reads))
  cleaned_reads$fastq_cleaned_reads <- as.numeric(as.character(cleaned_reads$fastq_cleaned_reads))
  read_stats <- data.frame (sample = colnames(bwa_counts_unique)[2:ncol(bwa_counts_unique)],
                            bwa_counts_unique = colSums(bwa_counts_unique[,2:ncol(bwa_counts_unique)]),
                            bwa_counts_unique_pass = colSums(bwa_counts_unique_filtered[,2:ncol(bwa_counts_unique_filtered)]),
                            bwa_counts_total = colSums(bwa_counts_total[,2:ncol(bwa_counts_total)]),
                            bwa_counts_total_pass = colSums(bwa_counts_total_filtered[,2:ncol(bwa_counts_total_filtered)]))
  read_stats <- merge(read_stats, raw_reads, by = "sample")
  read_stats <- merge(read_stats, cleaned_reads, by = "sample")
  read_stats$host_pct <- ((read_stats$fastq_raw_reads - read_stats$fastq_cleaned_reads) / read_stats$fastq_raw_reads) * 100
  read_stats$bwa_counts_unique_pct <- read_stats$bwa_counts_unique / read_stats$fastq_cleaned_reads * 100
  read_stats$bwa_counts_unique_pass_pct <- read_stats$bwa_counts_unique_pass / read_stats$fastq_cleaned_reads * 100
  read_stats$bwa_counts_total_pct <- read_stats$bwa_counts_total / read_stats$fastq_cleaned_reads * 100
  read_stats$bwa_counts_total_pass_pct <- read_stats$bwa_counts_total_pass / read_stats$fastq_cleaned_reads * 100
  
  # # Let's add in kneaddata map stats
  # if (cpu == "local") {
  #   metaphlan_reads <- read.csv("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Kneaddata_summary/Kneaddata_HG37_stats.csv")
  # } else if (cpu == "remote") {
  #   metaphlan_reads <- read.csv("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metamap/Kneaddata_HG37_stats.csv")
  # }
  # metaphlan_reads <- metaphlan_reads[!duplicated(metaphlan_reads$sampleid), ]
  # metaphlan_reads <- metaphlan_reads[,c(1,4:ncol(metaphlan_reads))]
  # colnames(metaphlan_reads) <- paste0("kneaddata_", colnames(metaphlan_reads))
  # colnames(metaphlan_reads)[1] <- "sample"
  # metaphlan_reads <- data.frame(lapply(metaphlan_reads, function(x) {gsub("%", "", x)}))
  # metaphlan_reads[, c(2:ncol(metaphlan_reads))] <- sapply(metaphlan_reads[, c(2:ncol(metaphlan_reads))], as.numeric)
  # metaphlan_reads <- metaphlan_reads[,c(1,2,9)]
  # read_stats <- merge(read_stats, metaphlan_reads, by = "sample")
  # colnames(read_stats)[14] <- "kneaddata_host_alignment"
  
  # Reorder dataframe
  #read_stats <- read_stats[, c((1:5),(9:12), (6:8), (13:14))]
  read_stats <- read_stats[, c((1:5),(9:12), (6:8))]
  
  # Write out tables
  if (cpu == "local") {
    write.table(species_richness_summary, paste0("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assembly comparison - species richness/speciesRichness_metaphlanVsmetamap_", prefix, "_", filter, "Filter.tsv"), sep = "\t", row.names = FALSE)
    write.table(read_stats, paste0("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assembly comparison - mapping rates/metamap_", prefix, "_", filter, "Filter_kneaddata_mapStats.tsv"), sep = "\t", row.names = FALSE)
  } else if (cpu == "remote") {
    write.table(species_richness_summary, paste0("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metamap/AssemblyComparison_speciesRichness/speciesRichness_metaphlanVsmetamap_", prefix, "_", filter, "Filter.tsv"), sep = "\t", row.names = FALSE)
    write.table(read_stats, paste0("/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metamap/AssemblyComparison_mappingRates/metamap_", prefix, "_", filter, "Filter_kneaddata_mapStats.tsv"), sep = "\t", row.names = FALSE)
  }
  
  # ## Let's plot species richness for metamap-UHGG vs metamap-metawrap
  # setwd("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap")
  # #species_richness_summary_UHGG <- read.table("speciesRichness_metaphlanVsmetamap_UHGG_strongFilter.tsv", header = TRUE)
  # species_richness_summary_UHGG <- read.table("speciesRichness_metaphlanVsmetamap_UHGG_weakFilter.tsv", header = TRUE)
  # #species_richness_summary_metawrap <- read.table("speciesRichness_metaphlanVsmetamapmetawrap.tsv", header = TRUE)
  # species_richness_summary_metawrap <- species_richness_summary
  # colnames(species_richness_summary_UHGG)[3:4] <- c("UHGG_total", "UHGG_unique")
  # colnames(species_richness_summary_metawrap)[3:4] <- c("metawrap_total", "metawrap_unique")
  # species_richness_summary_metawrap <- species_richness_summary_metawrap[,c(1,3,4)]
  # species_richness_summary <- merge(species_richness_summary_UHGG, species_richness_summary_metawrap, by = "sample")
  
  # # And plot
  # plot(species_richness_summary$UHGG_unique, species_richness_summary$metawrap_unique, pch=19)
  # plot(species_richness_summary$UHGG_total, species_richness_summary$metawrap_total, pch=19)
  # plot(species_richness_summary$metaphlan, species_richness_summary$metawrap_unique, pch=19)
  # plot(species_richness_summary$metaphlan, species_richness_summary$UHGG_unique, pch=19)
}