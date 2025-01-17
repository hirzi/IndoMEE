rm(list=ls())

# bio
# bio <- read.delim("input_files/clin_data_trim.tsv")
bio <- read.delim("input_files/Punan_160_full_ClinLs_metadata.tsv")

# fix gender
# bio[bio$Gender=="P","Gender"] <- "F"
# bio[bio$Gender=="L","Gender"] <- "M"

# obesity
bio$bmi <- round(bio$weight/((bio$height/100)^2),1)
bio$obesity <- ifelse(bio$bmi >30, "obese", "lean")

# assign population
# mypop <- read.table("input_files/MAGs/population_group.tsv", header=T)
# bio <- merge.data.frame(bio, mypop, by="sampleid")

# sample list
MAGs <- read.delim("input_files/MAGs/revised/Hybrid_wSPMP/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv", header = T)
sample.list <- colnames(MAGs)[27:263]
spore <- sample.list[grep("ERR", sample.list)]
neg.control <- sample.list[grep("HFC", sample.list)]
pos.control <- sample.list[grep("POS", sample.list)]
repeats <- sample.list[grep("B1|B2|B3|EREP", sample.list)]
dropped <- sample.list[grep("MALRPN001", sample.list)]

exclude <- unique(c(neg.control, pos.control, repeats, dropped, spore))
sampleid <- sample.list[!sample.list %in% exclude]


# filter table
df <- bio
df$population <- factor(df$population, levels = c("Asmat_DAI","Punan_HUL", "Punan_PBS", "Punan_SUL",
                                                  "Punan_RPN", "Basap_BRU", "Balinese_PDW", "Balinese_DPS"))

df <-  df[df$sampleid %in% sampleid, ]

write.table(df[,c("sampleid","population","Age","Gender","bmi","obesity", "Rice","Sago","Pork","Chicken","Fish","Egg","Hunting","Watching_TV")], 
            file="output_files/TableS1_2.Subject_Characteristics.tsv",
            row.names = F, quote = F, sep="\t")
write.table(df[,c("sampleid","population","Age","Gender","bmi","obesity", "Rice","Sago","Pork","Chicken","Fish","Egg","Hunting","Watching_TV")], 
            file = "output_files/MAGs/fig_out/Supplementary_Tables/TableS1_2.Subject_Characteristics_metadata.tsv", 
            row.names = F,quote = T,sep = "\t")


df <- df[,c("sampleid","population","Age","Gender","bmi", "obesity")]

apply(df, 1, function(x) any(is.na(x)))
apply(df[,c("Age","bmi")], 2, function(x) shapiro.test(x))

library(dplyr)

df %>% 
  select("sampleid","population","Age","Gender","bmi","obesity") %>%
  data.frame(stringsAsFactors = T) %>% 
  group_by(population) %>%
  summarise(n=length(sampleid),
            male=sum(sum(na.omit(Gender)=="M")),
            female=sum(sum(na.omit(Gender)=="F")),
            lean=sum(sum(na.omit(obesity)=="lean")),
            obese=sum(sum(na.omit(obesity)=="obese")),
            age=paste0(median(na.omit(Age))," (",min(na.omit(Age))," - ",max(na.omit(Age)),")"),
            BMI=paste0(median(na.omit(bmi))," (",min(na.omit(bmi))," - ",max(na.omit(bmi)),")"),
            Incomplete=paste0(ifelse(sum(is.na(Age))==0, "Age: ", paste0("Age:",sum(is.na(Age)))), " | ",
                              ifelse(sum(is.na(bmi))==0, "BMI: ", paste0("BMI:",sum(is.na(bmi)))), " | ",
                              ifelse(sum(is.na(Gender))==0, "Gender: ", paste0("Gender:",sum(is.na(Gender)))))) -> baseline_by_pop

df %>% 
  select("sampleid","population","Age","Gender","bmi","obesity") %>%
  data.frame(stringsAsFactors = T) %>% 
  summarise(population="All_Sample",
            n=length(sampleid),
            male=sum(sum(na.omit(Gender)=="M")),
            female=sum(sum(na.omit(Gender)=="F")),
            lean=sum(sum(na.omit(obesity)=="lean")),
            obese=sum(sum(na.omit(obesity)=="obese")),
            age=paste0(median(na.omit(Age))," (",min(na.omit(Age))," - ",max(na.omit(Age)),")"),
            BMI=paste0(median(na.omit(bmi))," (",min(na.omit(bmi))," - ",max(na.omit(bmi)),")"),
            Incomplete=paste0(ifelse(sum(is.na(Age))==0, "Age: ", paste0("Age:",sum(is.na(Age)))), " | ",
                              ifelse(sum(is.na(bmi))==0, "BMI: ", paste0("BMI:",sum(is.na(bmi)))), " | ",
                              ifelse(sum(is.na(Gender))==0, "Gender: ", paste0("Gender:",sum(is.na(Gender)))))) -> baseline


baseline_summary <- rbind(baseline_by_pop, baseline)

write.table(baseline_summary, file = "output_files/baseline_chara_summary.txt", row.names = F,quote = F,sep = "\t")
write.table(baseline_summary, file = "output_files/MAGs/fig_out/Supplementary_Tables/TableS1_1.Baseline_Summary.tsv", row.names = F,quote = T,sep = "\t")

saveRDS(sampleid, file = "output_files/sampleid_n116.rds")
saveRDS(spore, file = "output_files/sampleid_spore.rds")
exclude2 <- unique(c(neg.control, pos.control, repeats, dropped))
saveRDS(exclude2, file = "output_files/sampleid_exclusion.rds")
