rm(list=ls())

setwd("C:/Users/caf77/OneDrive - University of Cambridge/Documents/Analysis/MOBILE/0_Pilot_reviewed/")

library(vegan)
library(AICcPermanova)
library(ggplot2)
library(ggpubr)
library(ggtext)
source("fit_models_function.R")

##################
## Prepare Data ##
##################
# load metadata
env.df <- readRDS("output_files/revised_MAG/pmv.metadata.indospmp.rds")
env.df$cLifestyle <- ifelse(env.df$population %in% c("Asmat_DAI","Punan_PBS","Punan_HUL","Punan_SUL"), "remote",
                            ifelse(env.df$population %in% c("Punan_RPN","Basap_BRU","Balinese_PDW"), "rural",
                                   ifelse(env.df$population %in% c("Balinese_DPS","Chinese_SIN","Malay_SIN","Indian_SIN"),"urban", NA)))
env.df$lifestyle <- ifelse(env.df$population %in% c("Asmat_DAI","Punan_PBS","Punan_HUL","Punan_SUL"), 0,
                           ifelse(env.df$population %in% c("Punan_RPN","Basap_BRU","Balinese_PDW"), 0.5,
                                  ifelse(env.df$population %in% c("Balinese_DPS","Chinese_SIN","Malay_SIN","Indian_SIN"),1, NA)))

# age and gender
age_gender <- read.delim("input_files/input_age_sex_pop.tsv")[,c("sampleid","Age","Gender")]
env.df <- merge.data.frame(env.df, age_gender, by="sampleid")

# time since collection
sequencing_year <- 2021
year_collection <- read.delim("input_files/metadata_methodology.tsv")
year_collection$year_since_collection <- sequencing_year - year_collection$year_collection

# merge
env.df <- merge.data.frame(env.df, year_collection[,c("sampleid","year_since_collection","storage_method")], by="sampleid")
write.table(env.df, "output_files/revised_MAG/pmv.metadata.indospmp_revised.tsv", row.names = F, quote = F,sep = "\t")

# save metadata by dataset
# IndoMEE + SPMP
selectid <- rownames(readRDS("output_files/revised_MAG/pmv.comm.MAGs.indospmp.rds"))
select.metadata <- env.df
rownames(select.metadata) <- select.metadata$sampleid
select.metadata <- select.metadata[select.metadata$sampleid%in%selectid,]
select.metadata <- select.metadata[selectid,]
saveRDS(select.metadata,"output_files/revised_MAG/pmv.metadata.indospmp_revised.rds")


# IndoMEE Only
selectid <- rownames(readRDS("output_files/revised_MAG/pmv.comm.MAGs.indomee.rds"))
select.metadata <- env.df
rownames(select.metadata) <- select.metadata$sampleid
select.metadata <- select.metadata[select.metadata$sampleid%in%selectid,]
select.metadata <- select.metadata[selectid,]
saveRDS(select.metadata,"output_files/revised_MAG/pmv.metadata.indomee_revised.rds")


###############
## PERMANOVA ##
###############
# define ouput path
outdir="output_files/revised_MAG/PERMANOVA"

# model list
AllModels <- make_models(vars=c("lifestyle","Age","Gender","population","storage_method","year_since_collection"), ncores = 2)

# IndoMEE + SPMP
comm.df <- readRDS("output_files/revised_MAG/pmv.comm.MAGs.indospmp.rds")
env.df <- readRDS("output_files/revised_MAG/pmv.metadata.indospmp_revised.rds")
env.df$lifestyle <- as.numeric(env.df$lifestyle)
input.env <- env.df[,c("lifestyle","Age","Gender","population","storage_method","year_since_collection")]
input.env <- unique(data.frame(apply(env.df[,c("lifestyle","Age","year_since_collection")], 2,  function(x) as.numeric(x)),
                               env.df[,c("sampleid","Gender","population","storage_method")],
                               row.names = rownames(env.df)))
dim(input.env)

NonLinear <- filter_vif(all_forms = AllModels, env_data = input.env, filter = F, ncores = 3)
View(NonLinear)
write.table(NonLinear, file.path(outdir,"/PERMANOVA_CoLinearity_IndoSPMP.tsv"), quote = F, row.names = F, sep = "\t")

SPMP.Fit <- fit_models(all_forms = AllModels, 
                       veg_data = comm.df+1, 
                       env_data = input.env, 
                       method = "aitchison",
                       logfile = file.path(outdir,"/spmp_permanova.log"),
                       ncores = 3)

# rearrange & add collinearity
SPMP.Fit <- SPMP.Fit[order(SPMP.Fit$max_vif,decreasing = F),]
LS_select <- grep("lifestyle", SPMP.Fit$form)
SPMP.Fit$Model <- "no_lifestyle"
SPMP.Fit$Model[LS_select] <- "with_lifestyle"
SPMP.Fit$Model[SPMP.Fit$form == "Distance ~ 1"] <- "base"
SPMP.Fit <- rbind(SPMP.Fit[LS_select,], SPMP.Fit[-LS_select,])
SPMP.Fit$collinearity <- ifelse(SPMP.Fit$max_vif>=6, "high","low") # VIF cutoff at 6
write.table(SPMP.Fit, file.path(outdir,"/PERMANOVA_IndoSPMP.tsv"), quote = F, row.names = F, sep = "\t")

# aic_weight <- select_models(SPMP.Fit[SPMP.Fit$Model!="no_lifestyle",],delta_aicc = 10000)
aic_weight <- select_models(SPMP.Fit, delta_aicc = 1000)

print(aic_weight)
write.table(aic_weight, file.path(outdir,"/PERMANOVA_IndoSPMP_selected.tsv"), quote = F, row.names = F, sep = "\t")

masked <- which(colnames(aic_weight) %in% c("Model","collinearity"))
rsq.SPMP <- data.frame(akaike_adjusted_rsq(aic_weight[,-masked]))
print(rsq.SPMP)
write.table(rsq.SPMP, file.path(outdir,"/PERMANOVA_IndoSPMP_selected_AkaikeAdjR2.txt"), quote = F, row.names = F, sep = "\t")

# save input files for later
comm.spmp <- comm.df
env.spmp <- input.env
Model0 <- adonis2((comm.spmp+1) ~ lifestyle + storage_method, method="aitchison", env.spmp, by="margin", permutations = 9999, parallel = 3)
print(Model0)
saveRDS(Model, file.path(outdir,"/PERMANOVA_IndoSPMP_model_single_lifestyle+storage.rds"))

Model1 <- adonis2((comm.spmp+1) ~ lifestyle + Age + storage_method, method="aitchison", env.spmp, by="margin", permutations = 9999, parallel = 3)
print(Model1)
saveRDS(Model1, file.path(outdir,"/PERMANOVA_IndoSPMP_model_single_lifestyle+age+storage.rds"))


################
# IndoMEE Only #
################
comm.df <- readRDS("output_files/revised_MAG/pmv.comm.MAGs.indomee.rds")
env.df <- readRDS("output_files/revised_MAG/pmv.metadata.indomee_revised.rds")
env.df$lifestyle <- as.numeric(env.df$lifestyle)
input.env <- env.df[,c("lifestyle","Age","Gender","population","storage_method","year_since_collection")]
input.env <- unique(data.frame(apply(env.df[,c("lifestyle","Age","year_since_collection")], 2,  function(x) as.numeric(x)),
                               env.df[,c("sampleid","Gender","population","storage_method")],
                               row.names = rownames(env.df)))
dim(input.env)

NonLinear <- filter_vif(all_forms = AllModels, env_data = input.env, filter = F, ncores = 3)
View(NonLinear)
write.table(NonLinear, file.path(outdir,"/PERMANOVA_CoLinearity_IndoMEE.tsv"), quote = F, row.names = F, sep = "\t")


IndoMEE.Fit <- fit_models(all_forms = AllModels, 
                          veg_data = comm.df+1, 
                          env_data = input.env, 
                          method = "aitchison",
                          logfile = file.path(outdir,"/indomee_permanova.log"),
                          ncores = 2)

# rearrange & add collinearity
IndoMEE.Fit <- IndoMEE.Fit[order(IndoMEE.Fit$max_vif,decreasing = F),]
LS_select <- grep("lifestyle", IndoMEE.Fit$form)
IndoMEE.Fit$Model <- "no_lifestyle"
IndoMEE.Fit$Model[LS_select] <- "with_lifestyle"
IndoMEE.Fit$Model[IndoMEE.Fit$form == "Distance ~ 1"] <- "base"
IndoMEE.Fit <- rbind(IndoMEE.Fit[LS_select,], IndoMEE.Fit[-LS_select,])
View(IndoMEE.Fit)
IndoMEE.Fit$collinearity <- ifelse(IndoMEE.Fit$max_vif>=6, "high","low") # VIF cutoff at 6
write.table(IndoMEE.Fit, file.path(outdir,"/PERMANOVA_IndoMEE.tsv"), quote = F, row.names = F, sep = "\t")

# aic_weight <- select_models(IndoMEE.Fit[IndoMEE.Fit$Model!="no_lifestyle",],delta_aicc = 10000)
aic_weight <- select_models(IndoMEE.Fit,delta_aicc = 10000)
print(aic_weight)
write.table(aic_weight, file.path(outdir,"/PERMANOVA_IndoMEE_selected.tsv"), quote = F, row.names = F, sep = "\t")

masked <- which(colnames(aic_weight) %in% c("Model","collinearity"))
rsq.IndoMEE <- data.frame(akaike_adjusted_rsq(aic_weight[,-masked]))
print(rsq.IndoMEE)
write.table(rsq.IndoMEE, file.path(outdir,"/PERMANOVA_IndoMEE_selected_AkaikeAdjR2.txt"), quote = F, row.names = F, sep = "\t")

# single permanova with permutation
comm.indomee <- comm.df # save for later
env.indomee <- input.env # save for later
Model2 <- adonis2((comm.indomee+1) ~ lifestyle + storage_method, method="aitchison", env.indomee, by="margin", permutations = 9999, parallel = 3)
print(Model2)
saveRDS(Model2, file.path(outdir,"/PERMANOVA_IndoMEE_model_single_lifestyle+storage.rds"))

Model3 <- adonis2((comm.indomee+1) ~ lifestyle + Age + storage_method, method="aitchison", env.indomee, by="margin", permutations = 9999, parallel = 3)
print(Model3)
saveRDS(Model3, file.path(outdir,"/PERMANOVA_IndoMEE_model_single_lifestyle+age+storage.rds"))

#####################
## Summarise Table ##
#####################
library(tidyr)
library(broom)
library(dplyr)

var <- rsq.SPMP[!is.na(rsq.SPMP$Full_Akaike_Adjusted_RSq),"Variable"]
adj.rsq <- rsq.SPMP[!is.na(rsq.SPMP$Full_Akaike_Adjusted_RSq),"Full_Akaike_Adjusted_RSq"]
names(adj.rsq) <- var
perm <- t(tidy(Model1) %>% filter(term %in% var) %>% select(R2, p.value))
colnames(perm) <- var
out <- data.frame(value=c("adj.R2","R2","p.value"), rbind(adj.rsq,  perm), row.names = NULL)
out$dataset <- "IndoMEE + SPMP"
out$form <- "Distance ~ lifestyle + age + storage_method"

var <- rsq.IndoMEE[!is.na(rsq.IndoMEE$Full_Akaike_Adjusted_RSq),"Variable"]
adj.rsq <- rsq.IndoMEE[!is.na(rsq.IndoMEE$Full_Akaike_Adjusted_RSq),"Full_Akaike_Adjusted_RSq"]
names(adj.rsq) <- var
perm <- t(tidy(Model3) %>% filter(term %in% var) %>% select(R2, p.value))
colnames(perm) <- var
tmp <- data.frame(value=c("adj.R2","R2","p.value"), rbind(adj.rsq,  perm), row.names = NULL)
tmp$dataset <- "IndoMEE"
tmp$form <- "Distance ~ lifestyle + age + storage_method"

out <- rbind(out, tmp)
write.table(out, file.path(outdir,"/PERMANOVA_compare_dataset_lifestyle+age+storage.tsv"), quote = F, row.names = F, sep = "\t")


var <- c("lifestyle","storage_method")
adj.rsq <- rsq.SPMP[rsq.SPMP$Variable%in%var,"Full_Akaike_Adjusted_RSq"]
names(adj.rsq) <- var
perm <- t(tidy(Model0) %>% filter(term %in% var) %>% select(R2, p.value))
colnames(perm) <- var
out <- data.frame(value=c("adj.R2","R2","p.value"), rbind(adj.rsq,  perm), row.names = NULL)
out$dataset <- "IndoMEE + SPMP"
out$form <- "Distance ~ lifestyle + storage_method"

var <- rsq.IndoMEE[rsq.IndoMEE$Variable%in%var,"Variable"]
adj.rsq <- rsq.IndoMEE[!is.na(rsq.IndoMEE$Full_Akaike_Adjusted_RSq),"Full_Akaike_Adjusted_RSq"]
names(adj.rsq) <- var
perm <- t(tidy(Model2) %>% filter(term %in% var) %>% select(R2, p.value))
colnames(perm) <- var
tmp <- data.frame(value=c("adj.R2","R2","p.value"), rbind(adj.rsq,  perm), row.names = NULL)
tmp$dataset <- "IndoMEE"
tmp$form <- "Distance ~ lifestyle + storage_method"

out <- rbind(out, tmp)
write.table(out, file.path(outdir,"/PERMANOVA_compare_dataset_lifestyle+storage.tsv"), quote = F, row.names = F, sep = "\t")



