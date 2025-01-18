### This script converts metamap read counts to either RPKM and relative abundances.
## For reference, see: https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
# Method:
# 1. Normalise by total number of reads in sample (potentially divided by 1 millions for RPM). This normalises for sequencing depth.
# 2. Then, normalise by length of genome (by bp or by kb, the latter is known as RPKM)
# 3. To the get relative abundance, divide these values (RPKMs) by total sum of RPKMs.

# Set working directories
#counts_dir <- "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assemblies/On_hostDecontaminated_fastqs/"
counts_dir <- "/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metamap/metamap_on_metaWrap_MAGs_results/"
#metadata_dir <- "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assembly comparison - mapping rates/On_hostDecontaminated_fastqs/"
#metadata_dir <- "/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metamap/AssemblyComparison_mappingRates/On_hostDecontaminated_fastqs/"

# Set input (format) parameters
first_sample <- "ASMDAI005"

# Input assemblies
assembly_list <- list.dirs(counts_dir, recursive = FALSE, full.names= FALSE)
#assembly_list <- c("HybridAssembly_Comp50Cont5_draft2", "Hybrid_UHGG_SPMP_Assembly_Comp50Cont5_draft2", "SPMP_sANI95", "SPMP", "UHGG") #with summary_wSPMP
#prefix <- "HybridAssembly_Comp50Cont5_draft2"
#filt <- '*total*strongFilter.tsv*'

# List of filtering regimes
filter_list <- c('*total*strongFilter.tsv*', '*total*weakFilter.tsv*', '*unique*strongFilter.tsv*', '*unique*weakFilter.tsv*')

# Loop over assemblies
for (prefix in assembly_list) {
  
  for (i in seq(1,length(filter_list))) {
    
    # Define input (metamap counts output) files
    filt <- filter_list[i]
    summary_dir <- paste0(counts_dir, prefix, "/summary")
    #summary_dir <- paste0(counts_dir, prefix, "/summary_wSPMP")
    input <- list.files(summary_dir)[grepl(glob2rx(filt), list.files(summary_dir))]
    print(input)
    
    # Read in data
    counts_table_raw <- read.table(paste0(summary_dir, "/", input), header = TRUE)
    rownames(counts_table_raw) <- counts_table_raw$original_bin
    
    # Read in total reads per sample
    # if(grepl(glob2rx('*strong*'), filt) == TRUE) {
    #   mapstats_input <- list.files(metadata_dir)[grepl(glob2rx(paste0("*", prefix, "*strong*")), list.files(metadata_dir))]
    # } if(grepl(glob2rx('*weak*'), filt) == TRUE) {
    #   mapstats_input <- list.files(metadata_dir)[grepl(glob2rx(paste0("*", prefix, "*strong*")), list.files(metadata_dir))]
    # }
    # print(mapstats_input)
    # mapstats_raw <- read.table(mapstats_input, header = TRUE)
    # mapstats <- mapstats_raw[,c("sample", "fastq_raw_reads", "fastq_cleaned_reads")]
    
    # Format counts table
    colnames(counts_table_raw) <- gsub("X", "", colnames(counts_table_raw))
    sample_idx <- match(first_sample, colnames(counts_table_raw))
    if(is.na(match("length", colnames(counts_table_raw)))) {
      col_idx <- match("Length", colnames(counts_table_raw))
      colnames(counts_table_raw)[col_idx] <- "length"
    }
    length_idx <- match("length", colnames(counts_table_raw))
    #counts_table <- counts_table_raw[, c(1, length_idx, sample_idx:(ncol(counts_table_raw)-1))]
    counts_table <- counts_table_raw[, c(sample_idx:(ncol(counts_table_raw)-1))] # sample-specific
    #counts_table <- as.data.frame(counts_table_raw[, c(ncol(counts_table_raw))]) # sum across samples
    
    # Normalise by total number of reads in sample. For reference, see: https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
    counts_table_norm_reads <- apply(counts_table, 2, function(x) {x/sum(x)})
    
    # Normalise counts by length of genome, i.e. counts/length = counts/base pair = coverage
    counts_table_norm <- as.data.frame(sweep(as.matrix(counts_table_norm_reads), 1, counts_table_raw$length, `/`))
    
    # Convert to reads per kilobase million counts
    rpkm_counts_table <- counts_table_norm*1e9
    
    # Convert to relative abundance
    rel_abund_table <- as.data.frame(apply(counts_table_norm, 2, function(x) {x/sum(x)}))
    
    # Re-add and re-sort metadata to tables
    input_pref <- gsub(".tsv", "", input)
    rpkm_counts_table$SumAcrossAllSamples <- rowSums(rpkm_counts_table)
    if(all(rownames(counts_table_raw) == rownames(rpkm_counts_table)) == TRUE) {
      rpkm_counts_table <- cbind(counts_table_raw[,c(1:(sample_idx-1))], rpkm_counts_table)
      rpkm_counts_table <- rpkm_counts_table[order(rpkm_counts_table$SumAcrossAllSamples, decreasing = TRUE),]
      rpkm_counts_table <- rpkm_counts_table[,c(1:(ncol(rpkm_counts_table)-1))]
      rownames(rpkm_counts_table) <- seq(1,nrow(rpkm_counts_table))
      # Write out tables
      write.table(rpkm_counts_table, paste0(summary_dir,"/",input_pref,"_rpkm_revised.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
    } else {
      print("ERROR: Dataframes' orders don't match!")
    }
    rel_abund_table$SumAcrossAllSamples <- rowSums(rel_abund_table)
    if(all(rownames(counts_table_raw) == rownames(rel_abund_table)) == TRUE) {
      rel_abund_table <- cbind(counts_table_raw[,c(1:(sample_idx-1))], rel_abund_table)
      rel_abund_table <- rel_abund_table[order(rel_abund_table$SumAcrossAllSamples, decreasing = TRUE),]
      rel_abund_table <- rel_abund_table[,c(1:(ncol(rel_abund_table)-1))]
      rownames(rel_abund_table) <- seq(1,nrow(rel_abund_table))
      # Write out tables
      write.table(rel_abund_table, paste0(summary_dir,"/",input_pref,"_relAbund_revised.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
    } else {
      print("ERROR: Dataframes' orders don't match!")
    }
  }
}
