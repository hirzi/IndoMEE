# custom functions
convert_to_CPM <- function(comm.df){
  comm.df <- comm.df*10^6
  comm.df}
prevalence_filter <- function(comm.df, n) {
  x <- comm.df # matrix of abundance with sampleid as columns
  y <- x
  y[y>0] <- 1
  y <- y[rowSums(y)>n,]
  x <- x[rownames(x) %in% rownames(y),]
  x
}
headx <- function(data.frame,n=NULL){
  if(is.null(n)){
    data.frame[1:5,1:5]
  }else{
    data.frame[1:n,1:n]
  }
}

# load sampleid
sampleid <-readRDS("output_files/revised_MAG/pmv.sample.indospmp.rds")

# Metadata with Age and Gender
metadata <- readRDS("output_files/revised_MAG/pmv.metadata.indospmp.rds")[,c("sampleid","population","urbanisation")]
bio.indomee <- readRDS("output_files/clindata_imputed.rds")[,c("sampleid","Age","Gender")]
bio.spmp <- read.delim("input_files/MAGs/PRJEB49168_Illumina_WGS_samples_metadata.txt")[,c("run_accession","age","gender")]
colnames(bio.spmp) <- c("sampleid","Age","Gender")
colnames(bio.spmp) <- c("sampleid","Age","Gender")
bio <- rbind(bio.indomee,bio.spmp)
bio <- bio[bio$sampleid %in% sampleid, ]
metadata <- merge.data.frame(metadata, bio, by="sampleid")
rownames(metadata) <- metadata$sampleid
metadata <- metadata[sampleid,]
metadata$Age <- as.integer(metadata$Age)
write.table(metadata, "output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/input_age_sex_pop.tsv", row.names = F, quote = F, sep="\t")

obesity <- readRDS("output_files/clindata_imputed.rds")[,c("sampleid","bmi", "obesity")]
obesity <- obesity[!is.na(obesity$bmi),]
Ob <- obesity$sampleid[obesity$obesity=="obese"]
Ln <- obesity$sampleid[obesity$obesity=="lean"]

metadata <- metadata[!rownames(metadata) %in% Ob, ] # remove obese samples
metadata <- metadata[order(metadata$population),] ;  sampleid <- rownames(metadata ) # define sample order by population
metadata$lifestyle <- ifelse(metadata$urbanisation=="remote", 0, 
                             ifelse(metadata$urbanisation=="rural", 0.5,
                                    ifelse(metadata$urbanisation=="urban",1, NA)))
pop.ord <- levels(metadata$population)

# load and filter community (prevalence > 20)
df <- read.delim("input_files/MAGs/revised/Hybrid_wSPMP/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv", header = T)
reads <- df[,sampleid]
rownames(reads) <- df$original_bin
reads <- convert_to_CPM(comm.df = reads)
reads <- prevalence_filter(comm.df = reads,n=20)
preval <- reads ; preval[preval > 0] <- 1 ;table(rowSums(preval)>20) ; min(rowSums(preval)) # check prevalence

# load taxonomy
full.taxonomy <- read.csv("output_files/revised_MAG/MAASLIN2/indomee_full_taxonomy.txt", header=T)
full.taxonomy$feature <- full.taxonomy$original_bin

###########
# MAASLIN #
###########

# set input data
# a tab-delimited file with samples as rows and metadata as columns
input.comm <- data.frame(t(reads[,sampleid]))
input.metadata <- metadata[sampleid,]
table(rownames(input.comm) == rownames(input.metadata)) # check sample order

print(paste("Number of MAGs put into Maaslin:", ncol(input.comm)))
print(paste("Number of Samples put into Maaslin:", nrow(input.comm)))

write.table(data.frame(sampleid=rownames(input.comm),input.comm, row.names = NULL), 
            "output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/input.comm_p20.tsv", 
            quote = F, row.names = F, sep = "\t")
write.table(input.metadata,
            "output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/input.metadata_p20.tsv", 
            quote = F, row.names = F, sep = "\t")

###############
## Fit Model ##
###############
library("BiocManager")
# BiocManager::install("Maaslin2", force=TRUE)
library(Maaslin2)

######################################
# fit data discrete (medium as base) #
######################################
fit <- Maaslin2(input_data = input.comm,
                input_metadata = input.metadata,
                output="output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/cont_p20/",
                fixed_effects = c("lifestyle","Age","Gender"),
                normalization = "CLR",
                transform = "None",
                plot_scatter = F,
               save_models=T)


d0 <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/cont_p20/all_results.tsv", header=T) 
df <- d0[d0$qval<0.05,]
df <- df[df$metadata=="lifestyle",]

# count significant taxon
signif.feature <- length(unique(df$feature))
nfeature <- length(unique(d0$feature))
signif.pct <- round(100*signif.feature/nfeature,1)
varname <- paste0(unique(df$value), collapse = ",")
print(paste0(signif.pct,"% features (",signif.feature," of ", nfeature, ") were found significantly associated with ",varname))


################
# IndoMEE only #
################
indomee_sid <- readRDS("output_files/sampleid_n116.rds")
indomee_metadata <- read.delim("input_files/Punan_160_full_ClinLs_metadata_imputed.tsv")
indomee_metadata <- merge.data.frame(indomee_metadata, input.metadata[,c("sampleid","lifestyle")], by = "sampleid", all.y = T)
indomee_metadata <- indomee_metadata[indomee_metadata$sampleid %in% indomee_sid,]
rownames(indomee_metadata) <- indomee_metadata$sampleid
indomee_sid <- indomee_metadata$sampleid

# input
df <- read.delim("input_files/MAGs/revised/Hybrid_wSPMP/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv", header = T)
reads <- df[,indomee_sid]
rownames(reads) <- df$original_bin
reads <- convert_to_CPM(comm.df = reads)
preval <- reads ; preval[preval > 0] <- 1 ;table(rowSums(preval)>2) ; min(rowSums(preval)) # check prevalence
reads <- reads[which(rowSums(preval)>2),]

indomee_input <- data.frame(t(reads[,indomee_sid]))

fit <- Maaslin2(input_data = indomee_input,
                input_metadata = indomee_metadata,
                output="output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/cont_p20_indomee2/",
                #fixed_effects = c("lifestyle"),
                fixed_effects = c("lifestyle","Age","Gender"),
                normalization = "CLR",
                transform = "None",
                min_abundance = 2, 
                min_prevalence = 0,
                plot_scatter = F,
                save_models=T)



