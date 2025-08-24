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
df <- read.delim("input_files/abundance_tables_raw/original/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv")
metadata <- readRDS("input_files/pmv.metadata.indomee_revised.rds")
bio <- readRDS("input_files/clindata_imputed.rds")[,c("sampleid","Age","Gender","bmi","obesity")]
ffq <- read.delim("input_files/Punan_160_full_ClinLs_metadata.tsv",stringsAsFactors = T)

ls.var <- c("Rice","Sago","Tuber","Chicken","Pork","Fish","Egg","Milk","Hunting","Watching_TV")

###################
## Select Sample ##
###################
metadata <- metadata[metadata$sampleid %in% colnames(df),]
bio <- bio[bio$sampleid %in% colnames(df),c("sampleid","bmi","obesity")]
ffq <- ffq[ffq$sampleid %in% colnames(df),c("sampleid",ls.var)]

metadata <- merge.data.frame(metadata, bio, by="sampleid",all = T)
metadata <- merge.data.frame(metadata, ffq, by="sampleid", all = T)
metadata <- na.omit(metadata)

sampleid <- as.character(metadata$sampleid)

##################
## Filter Reads ##
##################
reads <- df[,sampleid]
rownames(reads) <- df$Genome
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
rownames(input.metadata) <- input.metadata$sampleid
input.metadata <- input.metadata[sampleid,]

###################
# Subset Taxonomy #
###################
source("0.parse_new_taxonomy_function.r")
full.taxonomy <- data.frame(df[,c("Genome","original_bin")], parse_taxonomy(df), stringsAsFactors = F) 
full.taxonomy$Species_deduplicated <- deduplicate_species(full.taxonomy$Species, type = "numeric")

subset.taxonomy <- full.taxonomy[full.taxonomy$Genome %in% colnames(input.comm),]

#write.table(full.taxonomy, "input_files/taxonomy/indomee_full_taxonomy.tsv", quote = F,sep = "\t",row.names = F)

################
## Save Input ##
################
write.table(input.comm, "output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_comm.tsv",quote = F,sep = "\t")
write.table(input.metadata,"output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_metadata.tsv",quote = F,sep = "\t")
write.table(full.taxonomy,"output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_full_taxonomy.tsv",quote = F,sep = "\t")

##################
## Reload Input ##
##################
# run new data
base_dir <- "C:/Users/caf77/OneDrive - University of Cambridge/Documents/Analysis/MOBILE/0_Pilot_reviewed/output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/"
input.comm <- read.delim(paste0(base_dir,"maaslin_input_comm.tsv", collapse = ""))
input.metadata <- read.delim(paste0(base_dir,"maaslin_input_metadata.tsv", collapse = ""))

input.metadata <- metadata[metadata$sampleid%in%rownames(input.comm),] 
rownames(input.metadata) <- input.metadata$sampleid
input.metadata <- input.metadata[rownames(input.comm),]

##################
# fit univariate #
##################
library(BiocManager)
library(Maaslin2)

output_dir0 <- "output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/"

dir.create(output_dir0,recursive = T,showWarnings = F)
adjustment <- c("Age", "Gender","bmi")
confounder <- c("year_since_collection", "storage_method")

# set Basap as the base factor for population
pop.ord <- as.character(unique(input.metadata$population))
n_basap <- which(pop.ord=="Basap_BRU")
pop.ord <- c(pop.ord[n_basap],pop.ord[-n_basap])
input.metadata$population <- factor(input.metadata$population, pop.ord)
levels(input.metadata$population )
print("Basap_BRU was used as the reference population.")

# for(i in c("population", ls.var)){ # switch on for the full run
for(i in ls.var[!ls.var %in% c("Sago","Pork","Chicken","population")]){ 
# for(i in c("Sago","Pork","Chicken","population")){ 
# for(i in c("Sago","Pork","Chicken")){ 
# for(i in "population"){  
  print(paste("Started MAASLIN2 for", i))
  
  output_dir_base <- paste0(output_dir0,i,"_noranef_base")
  output_dir_wpop_base <- paste0(output_dir0,i,"_popranef_base")
  
  fit <- Maaslin2(input_data = input.comm,
                  input_metadata = input.metadata,
                  output=output_dir_base,
                  fixed_effects = c(i),
                  normalization = "CLR",
                  transform = "None",
                  plot_scatter = F,
                  save_models=T,
                  cores = 4)
  
  if(i != "population"){
    fit <- Maaslin2(input_data = input.comm,
                    input_metadata = input.metadata,
                    output=output_dir_wpop_base,
                    fixed_effects = c(i),
                    random_effects = c("population"),
                    normalization = "CLR",
                    transform = "None",
                    plot_scatter = F,
                    save_models=T,
                    cores = 4)
  }
  
  output_dir <- paste0(output_dir0,i,"_noranef")
  output_dir_wpop <- paste0(output_dir0,i,"_popranef")

  fit <- Maaslin2(input_data = input.comm,
                  input_metadata = input.metadata,
                  output=output_dir,
                  fixed_effects = c(i, adjustment),
                  normalization = "CLR",
                  transform = "None",
                  plot_scatter = F,
                  save_models=T,
                  cores = 4)
  
  if(i != "population"){
    fit <- Maaslin2(input_data = input.comm,
                  input_metadata = input.metadata,
                  output=output_dir_wpop,
                  fixed_effects = c(i, adjustment),
                  random_effects = c("population"),
                  normalization = "CLR",
                  transform = "None",
                  plot_scatter = F,
                  save_models=T,
                  cores = 4)
  }
  
  
  output_dir_wConf <- paste0(output_dir0,i,"_noranef_wConf")
  output_dir_wConf_wpop <- paste0(output_dir0,i,"_popranef_wConf")
  
  fit <- Maaslin2(input_data = input.comm,
                  input_metadata = input.metadata,
                  output=output_dir_wConf,
                  fixed_effects = c(i, adjustment, confounder),
                  normalization = "CLR",
                  transform = "None",
                  plot_scatter = F,
                  save_models=T,
                  cores = 4)
  
  if(i != "population"){
  fit <- Maaslin2(input_data = input.comm,
                  input_metadata = input.metadata,
                  output=output_dir_wConf_wpop,
                  fixed_effects = c(i, adjustment, confounder),
                  random_effects = c("population"),
                  normalization = "CLR",
                  transform = "None",
                  plot_scatter = F,
                  save_models=T,
                  cores = 4)
  }
  
  print(paste(i, "finished running."))
}

#####################################
# Population with single-confounder #
#####################################
output_dir0 <- "output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/"

i="population"
output_dir_wConf2 <- paste0(output_dir0,i,"_noranef_wConf2")
fit <- Maaslin2(input_data = input.comm,
                input_metadata = input.metadata,
                output=output_dir_wConf2,
                fixed_effects = c(i, "storage_method"),
                normalization = "CLR",
                transform = "None",
                plot_scatter = F,
                save_models=T,
                cores = 4)

summary(fit$fits$MGYG000000921)

i="population"
output_dir_wConf3 <- paste0(output_dir0,i,"_noranef_wConf3")
fit2 <- Maaslin2(input_data = input.comm,
                input_metadata = input.metadata,
                output=output_dir_wConf3,
                fixed_effects = c(i, "year_since_collection"),
                normalization = "CLR",
                transform = "None",
                plot_scatter = F,
                save_models=T,
                cores = 4)

summary(fit2$fits$MGYG000000921)

#############################
# Population - NonBali Only #
#############################
output_dir0 <- "output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/"

nonBali_only <- as.character(input.metadata$sampleid[input.metadata$population %in% c("Asmat_DAI","Basap_BRU","Punan_HUL")])

# make sure we are filtering the same mags with the full model
prev_filter = 0.1*nrow(input.comm) # 0.1 is the maaslin default we used with the full Indomee dataset (n=73)
comm.df <- input.comm
comm.df[comm.df>0] <- 1
filtered_feature <- colSums(comm.df) > prev_filter
print(table(filtered_feature))

i="population"
output_dir3 <- paste0(output_dir0,i,"_noranef_nonBali_base")
fit3 <- Maaslin2(input_data = input.comm[nonBali_only,filtered_feature],
                 input_metadata = input.metadata[nonBali_only,],
                 output=output_dir3,
                 fixed_effects = c(i),
                 normalization = "CLR",
                 transform = "None",
                 plot_scatter = F,
                 save_models=T,
                 cores = 4,
                 min_prevalence = 0)

summary(fit3$fits$MGYG000000921) # check if model format is correct

i="population"
output_dir4a <- paste0(output_dir0,i,"_noranef_nonBali_adj")
fit4a <- Maaslin2(input_data = input.comm[nonBali_only,filtered_feature],
                  input_metadata = input.metadata[nonBali_only,],
                  output=output_dir4a,
                  fixed_effects = c(i, "Age","Gender","bmi"),
                  normalization = "CLR",
                  transform = "None",
                  plot_scatter = F,
                  save_models=T,
                  cores = 4,
                  min_prevalence = 0)

summary(fit4a$fits$MGYG000000921) # check if model format is correct

i="population"
output_dir4b <- paste0(output_dir0,i,"_noranef_nonBali_wConf")
fit4b <- Maaslin2(input_data = input.comm[nonBali_only,filtered_feature],
                  input_metadata = input.metadata[nonBali_only,],
                  output=output_dir4b,
                  fixed_effects = c(i, "Age","Gender","bmi","year_since_collection"),
                  normalization = "CLR",
                  transform = "None",
                  plot_scatter = F,
                  save_models=T,
                  cores = 4,
                  min_prevalence = 0)

summary(fit4b$fits$MGYG000000921) # check if model format is correct

# prelim mag count
signif_count <- data.frame()
signif_list <- list()

adjustment <- c("Age", "Gender", "bmi")
confounder <- c("year_since_collection", "storage_method")

#for(i in c("Sago", "Pork", "Chicken", "population")) {
for(i in c(ls.var, "population")) {
  output_dir0 <- "output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/"
  # output_dir0 <- "output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR_olddata/"
  
  if(i == "population") {
    suff <- c("_noranef_base", "_noranef", "_noranef_wConf", "_noranef_wConf2", "_noranef_wConf3")
  } else {
    suff <- c("_noranef_base", "_noranef", "_noranef_wConf", "_popranef_base", "_popranef", "_popranef_wConf")
  }
  
  for(m in suff){
    file_path <- paste0(output_dir0, i, m, "/all_results.tsv")
    
    if (!file.exists(file_path)) next  # Skip if file doesn't exist
    
    input.file <- read.delim(file_path)
    input.file <- input.file[input.file$qval < 0.05, ]
    input.file <- input.file[!input.file$metadata %in% c("Age", "Gender", "bmi", "year_since_collection", "storage_method"), ]
    
    ipop <- unique(input.file$metadata)
    if(all(ipop != "population")) {
      input.file <- input.file[input.file$metadata != "population", ]
    }
    
    signif.features <- unique(input.file$feature)
    
    if(length(signif.features) == 0){
      signif.features <- c("No Associations")
      out <- data.frame(model = paste0(i, m),
                        var = i,
                        ranef = gsub(pattern = "_", replacement = "", m),
                        signif.feature = 0)
    } else {
      out <- data.frame(model = paste0(i, m),
                        var = i,
                        ranef = gsub(pattern = "_", replacement = "", m),
                        signif.feature = length(signif.features))
    }
    
    out$ranef <- gsub("wConf", "_wConf", out$ranef)
    out$ranef <- gsub("base", "_base", out$ranef)
    
    # Define adjustment column
    if(m %in% c("_noranef_base", "_popranef_base")) {
      out$adjustment <- "none"
    } else if(m %in% c("_noranef", "_popranef")) {
      out$adjustment <- paste(adjustment, collapse = "+")
    } else if(m %in% c("_noranef_wConf", "_popranef_wConf")) {
      out$adjustment <- paste(c(adjustment, confounder), collapse = "+")
    } else if(m == "_noranef_wConf2") {
      out$adjustment <- "storage_method"
    } else if(m == "_noranef_wConf3") {
      out$adjustment <- "year_since_collection"
    }
    
    signif_count <- rbind(signif_count, out)
    
    nm <- paste0(i, m)
    signif_list[[nm]] <- signif.features
  }
}

# Display significant results only
signif_count[signif_count$signif.feature != 0 & signif_count$adjustment=="Age+Gender+bmi", ]

write.table(signif_count,"output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/0_SUMMARY/count_signif_all_forms.tsv",quote = F,sep = "\t",row.names = F)
# write.table(signif_count,"output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR_olddata/0_SUMMARY/count_signif_all_forms.tsv",quote = F,sep = "\t",row.names = F)

#################
# Multivariate ##
#################
base_dir <- "C:/Users/caf77/OneDrive - University of Cambridge/Documents/Analysis/MOBILE/0_Pilot_reviewed/output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/"

# set var
multivar <- c("Sago","Pork","Chicken") # we dropped Egg and Hunting because only 1 MAG was significant in the models with adjustments
popVar <- "population"
adjustmets <- c("Age", "Gender", "bmi")
confounders <- c("year_since_collection", "storage_method")
nonBali_only <- as.character(input.metadata$sampleid[input.metadata$population %in% c("Asmat_DAI","Basap_BRU","Punan_HUL")])
output_dir0 <- base_dir

for(sub in c("indomee","nonBali")){
  
  # prefilter features by prevalence
  prev_filter = 0.1*nrow(input.comm)
  comm.df <- input.comm
  comm.df[comm.df>0] <- 1
  filtered_feature <- colSums(comm.df) > prev_filter
  print(table(filtered_feature))
  
  if(sub == "indomee"){
    comm.df <- input.comm[,filtered_feature]
    meta.df <- input.metadata
    
    multi_dir0 <- paste0(output_dir0,"Multivar_diet")
    multi_dir1 <- paste0(output_dir0,"Multivar_diet+age+sex+bmi")
    multi_dir2 <- paste0(output_dir0,"Multivar_diet+age+sex+bmi+years+storage")
    
  } else if(sub == "nonBali"){
    comm.df <- input.comm[nonBali_only,filtered_feature]
    meta.df <- input.metadata[nonBali_only,]
    
    multi_dir0 <- paste0(output_dir0,"Multivar_diet_nonBali")
    multi_dir1 <- paste0(output_dir0,"Multivar_diet+age+sex+bmi_nonBali")
    multi_dir2 <- paste0(output_dir0,"Multivar_diet+age+sex+bmi+years_nonBali")
  }
  
  fit <- Maaslin2(input_data = comm.df,
                  input_metadata = meta.df,
                  output = multi_dir0,
                  fixed_effects = c(multivar),
                  normalization = "CLR",
                  transform = "None",
                  plot_scatter = F,
                  save_models=T,
                  min_prevalence = 0)
  
  fit <- Maaslin2(input_data = comm.df,
                  input_metadata = meta.df,
                  output = paste0(multi_dir0,"_popranef",collapse = ""),
                  fixed_effects = c(multivar),
                  random_effects = "population",
                  normalization = "CLR",
                  transform = "None",
                  plot_scatter = F,
                  save_models=T,
                  min_prevalence = 0)
  
  fit <- Maaslin2(input_data = comm.df,
                  input_metadata = meta.df,
                  output = multi_dir1,
                  fixed_effects = c(multivar, adjustmets),
                  normalization = "CLR",
                  transform = "None",
                  plot_scatter = F,
                  save_models=T,
                  min_prevalence = 0)
  
  fit <- Maaslin2(input_data = comm.df,
                  input_metadata = meta.df,
                  output = paste0(multi_dir1,"_popranef",collapse = ""),
                  fixed_effects = c(multivar, adjustmets),
                  random_effects = "population",
                  normalization = "CLR",
                  transform = "None",
                  plot_scatter = F,
                  save_models=T,
                  cores = 4, 
                  min_prevalence = 0)
  
  if(sub == "indomee"){
  fit <- Maaslin2(input_data = comm.df,
                  input_metadata = meta.df,
                  output = multi_dir2,
                  fixed_effects = c(multivar, adjustmets, "year_since_collection", "storage_method"),
                  normalization = "CLR",
                  transform = "None",
                  plot_scatter = F,
                  save_models=T,
                  min_prevalence = 0)
  } else {
    fit <- Maaslin2(input_data = comm.df,
                    input_metadata = meta.df,
                    output = multi_dir2,
                    fixed_effects = c(multivar, adjustmets, "year_since_collection"),
                    normalization = "CLR",
                    transform = "None",
                    plot_scatter = F,
                    save_models=T,
                    min_prevalence = 0)
  }
}

# population as fixed #
for(sub in c("indomee","nonBali")){
  
  # prefilter features by prevalence
  prev_filter = 0.1*nrow(input.comm)
  comm.df <- input.comm
  comm.df[comm.df>0] <- 1
  filtered_feature <- colSums(comm.df) > prev_filter
  print(table(filtered_feature))
  
  if(sub == "indomee"){
    comm.df <- input.comm[,filtered_feature]
    meta.df <- input.metadata
    
    multi_dir1 <- paste0(output_dir0,"Multivar_diet+age+sex+bmi+population")

  } else if(sub == "nonBali"){
    comm.df <- input.comm[nonBali_only,filtered_feature]
    meta.df <- input.metadata[nonBali_only,]
    
    multi_dir1 <- paste0(output_dir0,"Multivar_diet+age+sex+bmi+population_nonBali")
  }
  
  fit <- Maaslin2(input_data = comm.df,
                  input_metadata = meta.df,
                  output = multi_dir1,
                  fixed_effects = c(multivar, adjustmets),
                  normalization = "CLR",
                  transform = "None",
                  plot_scatter = F,
                  save_models=T,
                  min_prevalence = 0)
}


# count mags
base_dir <- "C:/Users/caf77/OneDrive - University of Cambridge/Documents/Analysis/MOBILE/0_Pilot_reviewed/output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/"
tests_list <- list.files(base_dir)
n <- c(grep("Multi",tests_list),
       grep("population",tests_list))
selected_tests <- unique(tests_list[n])

multi_list <- c(paste(base_dir, selected_tests,"/all_results.tsv",sep = ""),
                paste0(old_dir,"/Sago+Chicken+Pork_noranef","/all_results.tsv",collapse = ""),
                paste0(old_dir,"/population_noranef","/all_results.tsv", collapse = ""))

tmp <- lapply(multi_list, function(x) read.delim(x))
names(tmp) <- c(selected_tests,"old_base","old_population")

multi_list <- lapply(tmp, function(x) {
  x <- x[!x$metadata %in% c(adjustmets, confounder),]
  x <- x[x$qval < 0.05, ]
  return(x)
})  

multi_signif_list <- lapply(multi_list, function(x){
  x <- unique(x$feature)
  return(x)
})

multi_signif.cout <- data.frame()
for(i in 1:length(multi_list)){
  input.file <- multi_list[[i]]
  
  input.file <- input.file[input.file$qval < 0.05, ]
  input.file <- input.file[!input.file$metadata %in% c("Age", "Gender", "bmi", "year_since_collection", "storage_method"), ]
  signif.features <- unique(input.file$feature)
  
  m <- names(multi_list)[i]
  v <- unique(input.file$value)
  st <- length(grep("population", m)) != 0
  
  if(length(signif.features) == 0){
    signif.features <- c("No Associations")
    out <- data.frame(model = m,
                      var=ifelse(st==T, "population",paste0(c("Sago","Pork","Chicken"),collapse = "+")),
                      signif = paste0(v, collapse = ","),
                      signif.feature = 0)
  } else {
    out <- data.frame(model = m,
                      var=ifelse(st==T, "population",paste0(c("Sago","Pork","Chicken"),collapse = "+")),
                      signif = paste0(v, collapse = ","),
                      signif.feature = length(signif.features))
  }
  
  # Define adjustment column
  if(m %in% c("Multivar_diet", "Multivar_diet_nonBali", "population_noranef_base","population_noranef_nonBali_base","old_base", "old_population", "Multivar_diet_popranef","Multivar_diet_nonBali_popranef")) {
    out$adjustment <- "none"
  } else if (m %in% c("Multivar_diet+age+sex+bmi","Multivar_diet+age+sex+bmi","Multivar_diet+age+sex+bmi_nonBali","population_noranef","population_noranef_nonBali_adj","Multivar_diet+age+sex+bmi_popranef","Multivar_diet+age+sex+bmi_nonBali_popranef")){
    out$adjustment <- paste0(c("Age","Gender","BMI"), collapse = "+")
  } else if (m %in% c("Multivar_diet+age+sex+bmi+population","Multivar_diet+age+sex+bmi+population_nonBali")){
    out$adjustment <- paste0(c("Age","Gender","BMI","population"), collapse = "+") 
  } else if (m == "Multivar_diet+age+sex+bmi+years_nonBali") {
    out$adjustment <- paste0(c("Age","Gender","BMI","year_since_collection"), collapse = "+") 
  } else {
    out$adjustment <- paste0(c("Age","Gender","BMI","year_since_collection","storage_method"), collapse = "+")
  }
  
  multi_signif.cout <- rbind(multi_signif.cout, out)
  
}

multi_signif.cout$adjustment <- factor(multi_signif.cout$adjustment, c("none","Age+Gender+BMI","Age+Gender+BMI+population","Age+Gender+BMI+year_since_collection+storage_method","Age+Gender+BMI+year_since_collection"))
multi_signif.cout <- multi_signif.cout[order(multi_signif.cout$adjustment,decreasing = F),]
multi_signif.cout$ranef <- "none"
multi_signif.cout$ranef[grep("popranef",multi_signif.cout$model)] <- "popranef"
multi_signif.cout <- multi_signif.cout[order(multi_signif.cout$ranef),]
multi_signif.cout
write.table(multi_signif.cout, 
            file=paste0(output_dir0,"0_SUMMARY/count_signif_multivar_wnonBali.tsv", collapse = ""),
            quote = F, row.names = F, sep = "\t")

# plot
library(ggplot2)
library(ggvenn)
library(ggpubr)

venn_old <- multi_signif_list[c("Multivar_diet","population_noranef_base")]
venn_base <- multi_signif_list[c("Multivar_diet_nonBali","population_noranef_nonBali_base")]
venn_adjs <- multi_signif_list[c("Multivar_diet+age+sex+bmi","population_noranef")]
venn_conf <- multi_signif_list[c("Multivar_diet+age+sex+bmi_nonBali","population_noranef_nonBali_adj")]
venn_adjs_ranef <- multi_signif_list[c("Multivar_diet+age+sex+bmi_popranef","Multivar_diet+age+sex+bmi_nonBali_popranef","population_noranef")]

a <- ggvenn(venn_old, fill_color = c("tomato2","#7A67EE"), text_size = 7) + ggtitle("Base")
b <- ggvenn(venn_base, fill_color = c("tomato2","#7A67EE"), text_size = 7) + ggtitle("Base - NonBali")
c <- ggvenn(venn_adjs, fill_color = c("tomato2","#7A67EE"), text_size = 7) + ggtitle("Adjusted")
d <- ggvenn(venn_conf, fill_color = c("tomato2","#7A67EE"), text_size = 7) + ggtitle("Adjusted - NonBali")

ggarrange(plotlist = list(a,b,c,d))

