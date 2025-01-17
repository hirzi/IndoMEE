rm(list=ls())
setwd("C:/Users/caf77/OneDrive - University of Cambridge/Documents/Analysis/MOBILE/")

# setup
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(vegan)

pop.ord <- c("Asmat_DAI","Punan_PBS","Punan_HUL","Punan_SUL","Punan_RPN","Basap_BRU","Balinese_PDW","Balinese_DPS")
select.var <- c("Rice","Sago","Pork","Chicken","Fish","Egg","Milk","Hunting","Watching_TV")

mycols <- read.delim("input_files/MAGs/MAGs_pop_mycols.tsv", header=T)
mycols <- subset(mycols, group %in% pop.ord)
tmp <- mycols$mycols ; names(tmp) <- mycols$group
mycols <- tmp

dat <- read.table("input_files/final_metadata_ClinLs.tsv")
diet <- dat[,c("sampleid","population",select.var[1:7])]
tmp <- melt(diet, id.vars = c("sampleid","population"))
tmp <- na.omit(tmp)
tmp$population <- factor(tmp$population, pop.ord)


ggdata <- tmp
ggplot(ggdata, aes(x=variable, y=value, col=population, fill=population))+
  theme_pubclean(base_size = 20)+ theme(axis.title = element_blank())+
  geom_violin(scale = "width", alpha=0.5,trim = T,position = position_dodge(width = 0.8)) + scale_fill_manual(values = mycols) + scale_colour_manual(values = mycols) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8,jitter.width = 0.5,jitter.height = 0))+
  ggtitle("Annual consumption frequency") -> diet.freq

diet.freq

fileout="output_files/revised_MAG/FigureS3/FigS3A"
for(d in c("svg","pdf","png")){
  ggsave(plot=diet.freq, filename = paste0(fileout,".",d),device = d,width = 12,height = 8,units = "in")
}
saveRDS(diet.freq, file = paste0(fileout,".rds"))

actdata <- dat[,c("sampleid","population","Hunting","Watching_TV")]
actdata <- melt(actdata, id.vars = c("sampleid","population"),variable.name = "activity")
actdata <- na.omit(actdata) 
actdata %>% 
  group_by(activity,population) %>%
  summarise(yes=sum(value==1),
            no=sum(value==0)) %>%
  melt(id.vars = c("population","activity"),variable.name = "answer",value.name = "count") -> actdata

actdata$answer <- factor(actdata$answer,c("no","yes"))

ggplot(actdata)+ theme_pubclean(base_size = 14)+
  theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(legend.position = "none")+
  geom_bar(aes(x=population,y=count, fill=population, alpha=answer), col="grey30", 
           stat="identity", position = "stack")+
  scale_fill_manual(values = mycols) +
  scale_alpha_manual(values=c(0,1))+
  facet_wrap(.~activity,ncol = 1)+
  theme(axis.title.x = element_blank()) + ylab("frequency") -> act.plot
act.plot

fileout="output_files/revised_MAG/FigureS3/FigS3B"
for(d in c("svg","pdf","png")){
  ggsave(plot=act.plot, filename = paste0(fileout,".",d),device = d,width = 4,height = 8,units = "in")
}
saveRDS(act.plot, file = paste0(fileout,".rds"))
