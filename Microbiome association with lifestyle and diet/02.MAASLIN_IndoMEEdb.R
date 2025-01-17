rm(list=ls())

####################
# custom functions #
####################
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

# load data
df <- read.delim("input_files/MAGs/revised/Hybrid/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv")
metadata <- readRDS("output_files/revised_MAG/pmv.metadata.indomee.rds")[,c("sampleid","population","urbanisation")]
bio <- readRDS("output_files/clindata_imputed.rds")[,c("sampleid","Age","Gender","bmi","obesity")]
ffq <- read.delim("input_files/Punan_160_full_ClinLs_metadata.tsv",stringsAsFactors = T)
ls.var <- c("Rice","Sago","Tuber","Chicken","Pork","Fish","Egg","Milk","Hunting","Watching_TV")

###################
## Select Sample ##
###################
metadata <- metadata[metadata$sampleid %in% colnames(df),]
bio <- bio[bio$sampleid %in% colnames(df),c("sampleid","Age","Gender")]
ffq <- ffq[ffq$sampleid %in% colnames(df),c("sampleid",ls.var)]

metadata <- merge.data.frame(metadata, bio, by="sampleid")
metadata <- merge.data.frame(metadata, ffq, by="sampleid")

Ob <- bio$sampleid[bio$obesity=="obese"&!is.na(bio$obesity)]
incomplete <- as.character(ffq$sampleid)[apply(ffq,1, function(x) any(is.na(x)))]

sampleid <- as.character(metadata$sampleid)
sampleid <- sampleid[!sampleid %in% incomplete]
sampleid <- sampleid[!sampleid %in% Ob]

##################
## Filter Reads ##
##################
reads <- df[,sampleid]
rownames(reads) <- df$original_bin
reads <- convert_to_CPM(comm.df = reads)
reads <- prevalence_filter(comm.df = reads, n = 5) # abundance >0, prevalence >5
dim(reads)

print(paste("Number of complete cases:", ncol(reads))) # how many samples?
print(paste("Number of MAGs retained:", nrow(reads))) # how many MAGs?


########################
## Define Input Files ##
########################
# a tab-delimited file with samples as rows and metadata as columns

input.comm <- as.data.frame(t(reads))
input.comm <- input.comm[sampleid,]
input.metadata <- metadata[metadata$sampleid %in% sampleid, ]
rownames(input.metadata) <-input.metadata$sampleid
input.metadata <- input.metadata[sampleid,]

###############
## Fit Model ##
###############
library(Maaslin2)
print(paste("Using", paste(ls.var, collapse = ", "), "as fixed effect (ref. PERMANOVA Model Selection)."))

# set factor and level
# the first item is automatically set as reference
pop.ord <- as.character(unique(input.metadata$population))
n_basap <- which(pop.ord=="Basap_BRU")
pop.ord <- c(pop.ord[n_basap],pop.ord[-n_basap])
input.metadata$population <- factor(input.metadata$population, pop.ord)
print("Basap_BRU was used as the reference population.")

write.table(input.comm, "output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_comm.tsv",quote = F,sep = "\t")
write.table(input.metadata,"output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_metadata.tsv",quote = F,sep = "\t")

#########################
# fit univariate models #
#########################
output_dir0 <- "output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/"
ls.var

for(i in ls.var){
  output_dir <- paste0(output_dir0,i,"_noranef")
  output_dir_wpop <- paste0(output_dir0,i,"_popranef")
  
  print(paste("Started MAASLIN2 for", i))
  
  fit <- Maaslin2(input_data = input.comm,
                  input_metadata = input.metadata,
                  output=output_dir,
                  fixed_effects = c(i),
                  normalization = "CLR",
                  transform = "None",
                  plot_scatter = F,
                  save_models=T,
                  cores = 2)


  fit <- Maaslin2(input_data = input.comm,
                  input_metadata = input.metadata,
                  output=output_dir_wpop,
                  fixed_effects = c(i),
                  random_effects = c("population"),
                  normalization = "CLR",
                  transform = "None",
                  plot_scatter = F,
                  save_models=T,
                  cores = 2)
  
  print(paste(i, "finished running."))
}

# population as fixed effect
i="population"
output_dir <- paste0(output_dir0,i,"_noranef")

fit <- Maaslin2(input_data = input.comm,
                input_metadata = input.metadata,
                output=output_dir,
                fixed_effects = c(i),
                normalization = "CLR",
                transform = "None",
                plot_scatter = F,
                save_models=T,
                cores = 2)

# prelim mag count
signif_count <- data.frame()
signif_list <- list()
for(i in c(ls.var)){
  output_dir0 <- "output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/"
  for(m in c("_noranef","_popranef")){
    input.file <- paste0(output_dir0,i,m,"/all_results.tsv")
    input.file <- read.delim(input.file)
    input.file <- input.file[input.file$qval<0.05, ]
    signif.features <- unique(input.file$feature)
    
   
    if(length(signif.features)==0){
      signif.features <- c("No Associations")
      out <- data.frame(model=paste0(i,m), var=i, ranef=gsub(pattern = "_", replacement = "", m), signif.feature=0)
      signif_count <- rbind(signif_count, out)
    } else {
      out <- data.frame(model=paste0(i,m), var=i, ranef=gsub(pattern = "_", replacement = "", m), signif.feature=length(signif.features))
      signif_count <- rbind(signif_count, out)
    }
    
    len <- length(signif_list)
    
    if(len==0){
      signif_list[[1]] <- signif.features
      names(signif_list) <- paste0(i,m)
    } else {
      nm <- c(names(signif_list),paste0(i,m))
      signif_list[[len+1]] <- signif.features
      names(signif_list) <- nm
    }
  }
}

i="population"
output_dir0 <- "output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/"
m <- "_noranef"
input.file <- paste0(output_dir0,i,m,"/all_results.tsv")
input.file <- read.delim(input.file)
input.file$qval <- p.adjust(input.file$pval,method = "fdr")
input.file <- input.file[input.file$qval<0.05, ]
signif.features <- unique(input.file$feature)
if(length(signif.features)==0){
  signif.features <- c("No Associations")
  out <- data.frame(model=paste0(i,m), var=i, ranef=gsub(pattern = "_", replacement = "", m), signif.feature=0)
  signif_count <- rbind(signif_count, out)
  } else {
    out <- data.frame(model=paste0(i,m), var=i, ranef=gsub(pattern = "_", replacement = "", m), signif.feature=length(signif.features))
    signif_count <- rbind(signif_count, out)
    }
len <- length(signif_list)

nm <- c(names(signif_list),paste0(i,m))
signif_list[[len+1]] <- signif.features
names(signif_list) <- nm
      

signif_count[signif_count$signif.feature!=0,]

##############################
# select multivariate models #
##############################
multivar_list <- unique(signif_count[signif_count$signif.feature!=0,"var"])
exclude_var <- c("Hunting", "Watching_TV","population")
multivar_list <- multivar_list[!multivar_list%in% exclude_var]

AICcPermanova::make_models(multivar_list,k = 3) -> my.models

tmp <- unlist(lapply(strsplit(my.models$form ,split = "~"), function(x) x[2]))
model.list <- lapply(strsplit(tmp, split="\\+"), function(x) gsub(pattern = " ", replacement = "", x))

model.list[which(lapply(model.list, function(x) length(x)) >1)] -> model.list

output_dir0 <- "output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/"
for(i in 1:length(model.list)){
  ivar <- model.list[[i]]
  pref <- paste0(ivar,collapse = "+")
  output_dir <- paste0(output_dir0,pref,"_noranef")
  output_dir_wpop <- paste0(output_dir0,pref,"_popranef")
  
  print(paste("Started MAASLIN2 for", pref))
  
  fit <- Maaslin2(input_data = input.comm,
                  input_metadata = input.metadata,
                  output=output_dir,
                  fixed_effects = c(ivar),
                  normalization = "CLR",
                  transform = "None",
                  plot_scatter = F,
                  save_models=T,
                  cores = 2)
  
  
  fit <- Maaslin2(input_data = input.comm,
                  input_metadata = input.metadata,
                  output=output_dir_wpop,
                  fixed_effects = c(ivar),
                  random_effects = c("population"),
                  normalization = "CLR",
                  transform = "None",
                  plot_scatter = F,
                  save_models=T,
                  cores = 2)
  
  print(paste(pref, "finished running."))
}


for(i in 1:length(model.list)){
  ivar <- model.list[[i]]
  pref <- paste0(ivar,collapse = "+")
  
  output_dir0 <- "output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/"
  
  for(m in c("_noranef","_popranef")){
    input.file <- paste0(output_dir0,pref,m,"/all_results.tsv")
    input.file <- read.delim(input.file)
    input.file <- input.file[input.file$qval<0.05, ]
    signif.features <- unique(input.file$feature)

    if(length(signif.features)==0){
      signif.features <- c("No Associations")
      out <- data.frame(model=paste0(pref,m), var=pref, ranef=gsub(pattern = "_", replacement = "", m), signif.feature=0)
      signif_count <- rbind(signif_count, out)
    } else {
      out <- data.frame(model=paste0(pref,m), var=pref, ranef=gsub(pattern = "_", replacement = "", m), signif.feature=length(signif.features))
      signif_count <- rbind(signif_count, out)
    }
    
    len <- length(signif_list)
    nm <- c(names(signif_list),paste0(pref,m))
    signif_list[[len+1]] <- signif.features
    names(signif_list) <- nm
  }
  }



str(signif_list)
signif_count[signif_count$signif.feature>0,]

