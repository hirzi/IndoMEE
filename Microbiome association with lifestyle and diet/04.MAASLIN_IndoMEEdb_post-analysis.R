rm(list=ls())

# load input 
input.comm <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_comm.tsv")
dim(input.comm)
nMAG <- ncol(input.comm)
nSample <- nrow(input.comm)

# load taxonomy
full.taxonomy <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_full_taxonomy.tsv")

# Compare counts
significant.count <- readRDS("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/post_analysis/significant.MAGs.count.rds")

library(ggplot2)
library(ggpubr)
library(ggvenn)

test.ord <- unique(significant.count$var)
test.ord <- c(test.ord[grep("population",test.ord)],test.ord[-grep("population",test.ord)])
significant.count$var <- factor(significant.count$var, levels = test.ord)
significant.count$ranef <- factor(significant.count$ranef, c("noranef","popranef"))
# significant.count[significant.count$signif.feature!=0,]

ggdata <- significant.count[-grep(pattern = "Egg|Hunting|Watching_TV",significant.count$model),]

countbar <- ggplot(ggdata) + 
  theme_pubclean(base_size = 16) +
  geom_bar(aes(x=var, y=signif.feature, fill=ranef), stat = "identity", position = "dodge", col="black", width = 0.8) +
  geom_text(aes(x=var, y=signif.feature+10, col=ranef, label=signif.feature), position=position_dodge(width = 0.8))+
  scale_fill_manual(values=c("white","#7A67EE")) +
  scale_colour_manual(values=c("black","#7A67EE")) +
  scale_y_continuous(expand = c(0.01,0.01))+
  theme(axis.text.x=element_text(angle=90,vjust=0.5, hjust=1))+
  theme(axis.title=element_blank())+
  labs(title="Count of Significant MAGs", subtitle="Compare models with / without population as random effect")

countbar
outpath <- "output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/fig_out/"
for(d in c("svg","png","pdf")){
ggsave(plot = countbar, 
       filename = paste0(outpath,"signif_count_bar.",d), scale = 1.5,
       device = d,width = 10,height = 6,units = "in",dpi = "print")}


## Univariate Intersect ##
MAG.list <- readRDS("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/post_analysis/significant.MAGs.list.rds")

empties <- significant.count[significant.count$signif.feature==0,"model"]
model.names <- names(MAG.list)
model.names <- model.names[!model.names %in% empties]

noranef.list <- model.names[grep("noranef", model.names)]
noranef.list <- noranef.list[-grep("Hunting|Watching|Egg|Milk|Tuber|population", noranef.list)]

noranef_compare <- data.frame()
for(f in noranef.list){
  int <- intersect(x=MAG.list[[f]],
                   y=MAG.list[["population_noranef"]])
  AB <- length(int)
  A1 <- length(MAG.list[[f]]) - AB 
  B1 <- length(MAG.list[["population_noranef"]]) - AB
  AB1 <- nMAG - ( A1 + B1 + AB)
  
  out <- data.frame(model=f,
                    diet.only=A1,
                    intersect=AB,
                    pop.only=B1,
                    not.significant=AB1)
  noranef_compare <- rbind(noranef_compare, out)
}

noranef_compare_pct <- data.frame(model=noranef_compare$model,
                                  t(apply(noranef_compare[,-1],1, function(x) paste0("(",round(100*x/sum(x),1),"%)"))))
colnames(noranef_compare_pct) <- colnames(noranef_compare)

library(reshape2)
noranef_compare <- reshape2::melt(noranef_compare,id.vars = "model", variable.name = "set",value.name = "count")
noranef_compare_pct <- reshape2::melt(noranef_compare_pct,id.vars = "model", variable.name = "set",value.name = "label")
noranef_compare$label <- paste(noranef_compare$count, noranef_compare_pct$label,sep = "\n")

noranef_compare$model <- gsub(pattern = "_noranef",replacement = "", noranef_compare$model)
noranef_compare$model <- factor(noranef_compare$model, levels = test.ord)

noranef_compare$set <- factor(noranef_compare$set, rev(c("diet.only","intersect","pop.only","not.significant")))

noranef_compare <- noranef_compare[order(noranef_compare$set),]
noranef_compare <- noranef_compare[order(noranef_compare$model),]
noranef_compare$label <- factor(noranef_compare$label, unique(noranef_compare$label))

rownames(noranef_compare) <- NULL

popinter <- ggplot(noranef_compare) +
  theme_pubclean(base_size = 16) + 
  geom_bar(aes(x=model, y=count, fill=set), stat="identity", position="stack", width = 0.8) +
  geom_text(aes(x=model, y=count, label=label, group=set), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("grey90","grey80","#7A67EE","#FA8072")) +
  scale_y_continuous(expand = c(0.1,0.1))+
  theme(axis.title=element_blank()) +
  theme(axis.text.x=element_text(angle=90,vjust=0.5, hjust=1))+
  labs(title="Count of Significant MAGs", caption="Note: All models with no random effect. Population as fixed effect.")

popinter

outpath <- "output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/fig_out/"
for(d in c("svg","png","pdf")){
  ggsave(plot = popinter, 
         filename = paste0(outpath,"Population_intersect_count.",d), scale = 1.5,
         device = d,width = 10,height = 6,units = "in",dpi = "print")}

  
## Venn Analysis ###
diet.cols <- c("#2E8B57","#DB7093","#63B8FF")
#names(diet.cols) <- c("Sago","Pork","Chicken")

multivar.df <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/Sago+Chicken+Pork_noranef/significant_results_wTax.tsv")

vs <- list(Sago=multivar.df$feature[multivar.df$value=="Sago"],
           Pork=multivar.df$feature[multivar.df$value=="Pork"],
           Chicken=multivar.df$feature[multivar.df$value=="Chicken"])

p <- ggvenn(vs,text_size = 4, set_name_size = 5, fill_color = diet.cols)
print(p)

ggsave(p, file="output_files/revised_MAG/Figure3/Fig3C.svg",device = "svg",width = 8,height = 8)
ggsave(p, file="output_files/revised_MAG/Figure3/Fig3C.pdf",device = "pdf",width = 8,height = 8)
saveRDS(p, file="output_files/revised_MAG/Figure3/Fig3C.rds")



vs2 <- list(Population=MAG.list[["population_noranef"]],
            Diet=MAG.list[["Sago+Chicken+Pork_noranef"]])

p2 <- ggvenn(vs2,text_size = 4,set_name_size = 5, fill_color = c("#7A67EE","#FA8072"))
print(p2)

ggsave(p2, file="output_files/revised_MAG/Figure3/Fig3D.svg",device = "svg",width = 8,height = 8)
ggsave(p2, file="output_files/revised_MAG/Figure3/Fig3D.pdf",device = "pdf",width = 8,height = 8)
saveRDS(p2, file="output_files/revised_MAG/Figure3/Fig3D.rds")


##################
## diet volcano ##
##################
df <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/Sago+Chicken+Pork_noranef/all_results_wTax.tsv")
length(unique(df$feature)) -> ct
df <- subset(df, qval < 0.05)
length(unique(df$feature)) -> sig.ct

sig.ct *100 / ct

ggdata <- df[,c("feature","value","coef","N.not.0","pval","qval","species","Genome")]
ggdata$association <- ggdata$value
ggdata$association[ggdata$qval>=0.05] <- "none"
ggdata$association <- factor(ggdata$association, c("Pork","Sago","Chicken","none"))

ggdata <- ggdata[order(abs(ggdata$coef),decreasing = T),]
# ggdata$label <- seq(1,nrow(ggdata),by=1)
ggdata$label <- ggdata$species
ggdata$label[ggdata$qval>=0.05] <- ""
rownames(ggdata) <- NULL


diet.label <- ggdata[ggdata$label!="",]
diet.label <- merge.data.frame(diet.label, full.taxonomy[,c("feature","phylum","family","species")], by="feature")
rownames(diet.label) <- NULL
write.table(diet.label,"all_results(Sago+Pork+Chicken)_diet_label.tsv",
            quote = F, row.names = F, sep = "\t")

volc <- ggplot(ggdata, aes(y=-log10(qval), x=coef, fill=association)) + theme_bw() +
  geom_text(aes(label=label), col="grey70", size=3,nudge_y =0.1) +
  annotate("text", x = 0, y = -log10(0.045), label = "qval=0.05", vjust=0, col="red")+
  geom_hline(yintercept = -log10(0.05), col="red", lty=2)+
  geom_point(size=5, pch=21, col="grey30") +
  scale_fill_manual(values = c("#2E8B57","#FF69B4","#63B8FF","grey70")) +
  scale_x_continuous(limits = c(-3,3), breaks = sort(c(seq(-4,4, by=1),0)))+
  labs(title = "Diet Correlations (Multi-variate)",
       caption = "q-val is p-val adjusted with Benjamin Hochberg method, q-val<0.05 considered as significant") +
  ylab("Significance -log10(qval)") + xlab("Coefficient")

volc

outpath <- "output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/fig_out/"
for(d in c("svg","png","pdf")){
  ggsave(plot = volc, 
         filename = paste0(outpath,"diet_volcano.",d), scale = 1.5,
         device = d,width = 10,height = 6,units = "in",dpi = "print")}



#######################
## MAASLIN pop ranef ##
#######################
library(lme4)
library(reshape2)
library(dplyr)
library(ggplot2)
library(gridExtra)

# load model 
models <- readRDS("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/Sago+Chicken+Pork_popranef/fits/models.rds")
                  
pop.varex <- data.frame(feature=c(), pop_var=c())
for (i in names(models)){
  test <- models[[i]]
  ranef.var <- data.frame(VarCorr(test))[,c("grp","var1","vcov","sdcor")]
  varex <- ranef.var$vcov
  names(varex) <- ranef.var$grp
  varex <- round(varex*100/sum(varex),1)
  
  tmp <- data.frame(feature=i,pop_var=varex[1],row.names = NULL)
  
  pop.varex <- rbind(pop.varex,tmp)
}

# most abundant in
df.abund <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_comm.tsv",header = T)
df.abund$sampleid <- rownames(df.abund)
df.abund <- melt(df.abund, id.vars = "sampleid",variable.name = "feature",value.name = "CPM")

input.metadata <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_metadata.tsv",header = T)
input.metadata$sampleid <- rownames(input.metadata)
input.metadata <- input.metadata[,c("sampleid","population")]

df.abund <- merge.data.frame(input.metadata,df.abund, by="sampleid")
df.abund <- merge.data.frame(df.abund, full.taxonomy, by="feature")


df.abund %>% group_by(phylum, feature, population) %>% summarise(popMean=mean(CPM)) -> pop_average

mostabund.pop <- data.frame(feature=c(), topPop=c(), meanCPM=c())
for(i in unique(pop_average$feature)){
  tmp <- pop_average[pop_average$feature==i,]
  max <- max(tmp$popMean)
  n <- grep(max, tmp$popMean)
  
  j <- tmp$population[n]
  
  tmp <- data.frame(feature=i, topPop=j, meanCPM=max)
  mostabund.pop <- rbind(mostabund.pop, tmp)
}
mostabund.pop

popranef.df <- merge.data.frame(pop.varex, mostabund.pop, by="feature")
write.table(popranef.df, 
            "output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/Sago+Chicken+Pork_popranef/all_results(Sago+Pork+Chicken)_popvar.tsv", 
            quote = F,sep = "\t", row.names = F)

popvar.cutoff <- 35 # in percentage
strong.popranef <- unique(popranef.df[popranef.df$pop_var>=popvar.cutoff,"feature"])
length(strong.popranef)

###############
## diet list ##
###############
df <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/Sago+Chicken+Pork_popranef/all_results_wTax.tsv", header=T)

ggdata <- df[,c("feature","value","coef","N.not.0","pval","qval","phylum","class","family","species","Genome")]
ggdata <- ggdata[ggdata$qval<0.05,]
diet.MAGs <- unique(ggdata$feature)
length(diet.MAGs)

# order data
ggdata <- ggdata[order(ggdata$coef,decreasing = F),]
ggdata$value <- factor(ggdata$value, c("Pork","Sago","Chicken","none"))

mag.levels <- as.character(unique(ggdata$feature))
ggdata$feature <- factor(ggdata$feature, mag.levels)

phylum.tree.order <- rev(readRDS("output_files/revised_MAG/phylum.tree.order.rds"))
phylum.tree.order <- c(phylum.tree.order,"Methanobacteriota")
ggdata$phylum <- factor(ggdata$phylum, phylum.tree.order)

diet.list.df <- ggdata 

########
# plot #
########

#ggplot(diet.list.df, aes(x=coef,y=feature))+ theme_minimal() +
 # theme(plot.margin = margin(0,0,0,0),
  #      panel.spacing = unit(0, "cm"),
   #     legend.position = "right",
    #    panel.background = element_rect(color = "black",fill=NA),
     #   strip.text.y.left = element_text(angle=0),
      #  strip.background = element_rect(fill="white",linewidth = 0),
 #       text=element_text(size=12), axis.text=element_text(size=12),
  #      axis.text.y = element_text(angle=0, hjust = 1, size=8),
   #     # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    #    axis.title.y= element_blank(),
     #   panel.grid.major.x = element_line(colour="grey",linetype = 2)) + 
#  scale_color_manual(values = c("#00A600","#FF69B4","#0000FF","grey70")) +
 # scale_x_continuous(breaks = seq(-3,4, by=0.5),limits = c(-3,4), expand = c(0.05,0))+
  #xlab("Coefficients") +
#  labs(title="MAGs with Significant Diet Associations (qval<0.05)") +
 # facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
  #geom_point(aes(col=value),size=3) +
#  geom_text(aes(label=species), angle=0, nudge_x = 0.1, 
 #           size=3, col="grey40", hjust=0, vjust=0.2)

#
# abundance heatmap (total)
library(reshape2)
df.abund <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_comm.tsv",header = T)
df.abund$sampleid <- rownames(df.abund)
df.abund <- melt(df.abund, id.vars = "sampleid",variable.name = "feature",value.name = "CPM")

input.metadata <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_metadata.tsv",header = T)
input.metadata$sampleid <- rownames(input.metadata)
input.metadata <- input.metadata[,c("sampleid","population")]

df.abund <- merge.data.frame(input.metadata,df.abund, by="sampleid")
df.abund <- merge.data.frame(df.abund, full.taxonomy, by="feature")

n <- length(unique(df.abund$sampleid))

df.abund %>% group_by(feature, phylum) %>% summarise(abundance=sum(CPM)/(n*10^6)) -> total_abund

heat.df.total <- total_abund[total_abund$feature %in% mag.levels,]
heat.df.total$feature <- factor(heat.df.total$feature, mag.levels)
heat.df.total$phylum <- factor(heat.df.total$phylum, phylum.tree.order)


#ggplot(heat.df.total)+ theme_minimal() + 
#  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
#  geom_tile(aes(x="A", y=feature,fill=log(abundance)))+
#  scale_fill_gradient(low = "grey90", high = "black",)+
#  theme(plot.margin = margin(0,0,0,0),
#        panel.spacing = unit(0, "cm"),
 #       legend.position = "right",
  #      panel.grid = element_blank(),
   #     panel.border = element_rect(colour = "black", fill=NA),
    #    strip.text.y.left = element_text(angle=0),
     #   strip.background = element_rect(fill="white",linewidth = 0),
      #  text=element_text(size=12), axis.text=element_text(size=12),
       # axis.ticks.x= element_blank(),
#        axis.text.x=element_text(colour = "white"),
 #       # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
  #      axis.text.y = element_text(size=8),
   #     axis.title.y= element_blank()) +
#  labs(title="Total Abundance") +
#  scale_x_discrete(expand = c(0,0))+
#  scale_y_discrete(expand = c(0,0))+
#  xlab("")

# abundance heatmap (bypop)
# normalised average
# the normalisation is the percentage of the population average over the sum of the average value in each population by features
df.abund <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_comm.tsv",header = T)
df.abund$sampleid <- rownames(df.abund)
df.abund <- melt(df.abund, id.vars = "sampleid",variable.name = "feature",value.name = "CPM")

input.metadata <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_metadata.tsv",header = T)
input.metadata$sampleid <- rownames(input.metadata)
input.metadata <- input.metadata[,c("sampleid","population")]

df.abund <- merge.data.frame(input.metadata,df.abund, by="sampleid")
df.abund <- merge.data.frame(df.abund, full.taxonomy, by="feature")

# n <- length(unique(df.abund$sampleid))

df.abund %>% group_by(population,feature, phylum) %>% 
  summarise(n=length(unique(sampleid)),
            popCPM=sum(CPM),
            popMedian=median(CPM),
            popMean=mean(CPM),
            popSD=sd(CPM)) -> pop_abund

df.abund %>% group_by(feature) %>% summarise(totalCPM=sum(CPM)) -> total_abund

heat.df.pop <- merge.data.frame(pop_abund, total_abund, by="feature")

pop_abund %>% group_by(feature) %>% summarise(sum.popMean=sum(popMean)) -> sumpopMean
heat.df.pop <- merge.data.frame(heat.df.pop, sumpopMean, by="feature")

heat.df.pop <- heat.df.pop[order(heat.df.pop$feature),]
heat.df.pop$normalised.abund <- 100*heat.df.pop$popCPM/heat.df.pop$totalCPM
heat.df.pop$normalised.average <- 100*heat.df.pop$popMean/heat.df.pop$sum.popMean

# order data
heat.df.pop <- heat.df.pop[heat.df.pop$feature %in% mag.levels,]
heat.df.pop$feature <- factor(heat.df.pop$feature,mag.levels)
heat.df.pop$phylum <- factor(heat.df.pop$phylum, phylum.tree.order)
heat.df.pop$population <- factor(heat.df.pop$population,
                                 c("Asmat_DAI","Punan_HUL","Basap_BRU",
                                   "Balinese_PDW","Balinese_DPS"))

# ggplot(heat.df.pop, aes(x=population, y=feature, fill= normalised.average)) + 
 # facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
  #theme_minimal() + 
  #theme(plot.margin = margin(0,0,0,0),
#        panel.spacing = unit(0, "cm"),
#        legend.position = "right",
#        panel.grid = element_blank(),
#        panel.border = element_rect(colour = "black", fill=NA),
#        text=element_text(size=12), axis.text=element_text(size=12),
#        axis.ticks.x= element_blank(),
#        axis.text.x=element_text(colour = "black", size=8),
#        strip.text.y.left = element_text(angle=0),
#        strip.background = element_rect(fill="white",linewidth = 0),
#        # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#        axis.text.y = element_text(size=8),
#        axis.title.y= element_blank()) +
#  scale_x_discrete(expand = c(0,0))+
#  scale_y_discrete(expand = c(0,0))+
#  geom_tile()+ scale_fill_gradient(low="white",high="black") +
#  labs(title="Normalised Average") +
#  xlab("")

########
# plot #
########

diet.cols <- c("#00A600","#FF69B4","#0000FF","#7A67EE","grey70")
names(diet.cols) <- c("Sago","Pork","Chicken","population","not significant")

a <- ggplot(diet.list.df, aes(x=coef,y=feature))+ 
  theme_minimal() + 
  theme(plot.margin = margin(0,0,21,0,unit = "mm"),
        plot.background = element_blank(),
        panel.spacing = unit(0, "cm"),
        panel.grid=element_blank(),
        panel.grid.major.x = element_line(colour="grey",linetype = 3, linewidth = .5),
        panel.grid.minor.x = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.y.left = element_text(angle=0),
        strip.background = element_rect(fill="white",linewidth = 0),
        text=element_text(size=12), axis.text=element_text(size=12),
        #axis.text.y = element_text(angle=0, hjust = 1, size=8),
        axis.text.x=element_text(colour = "black", size=8),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y= element_blank())+ 
  scale_color_manual(values = diet.cols) +
  scale_x_continuous(breaks = seq(-3,4, by=1),limits = c(-3,4), expand = c(0.06,0.06))+
  scale_y_discrete(expand = c(0,0.5))+
  xlab("Coefficients") +
  labs(title="MAGs with Significant Diet Associations (qval<0.05)") +
  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
  geom_point(aes(col=value),size=3) +
  geom_text(aes(label=species), angle=0, nudge_x = 0.1, 
            size=3, col="grey40", hjust=0, vjust=0.2)

b <- ggplot(heat.df.total)+ 
  theme_minimal() +
  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
  geom_tile(aes(x="Total Abundance", y=feature,fill=log(abundance)))+
  scale_fill_gradient(low = "grey90", high = "black",)+
  theme(plot.margin = margin(0,0,0,0,unit = "mm"),
        plot.background = element_blank(),
        panel.spacing = unit(0, "cm"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        #strip.text.y.left = element_text(angle=0),
        #strip.background = element_rect(fill="white",linewidth = 0),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(size=10, colour = "black",angle=60, hjust = 1),
        axis.ticks.x= element_line(colour = NA),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        # axis.text.y = element_text(size=8),
        axis.title.y= element_blank()) +
  labs(title="") +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0.5))+
  xlab("")

c <- ggplot(heat.df.pop, aes(x=population, y=feature, fill=normalised.average)) + 
  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
  theme_minimal() + 
  theme(plot.margin = margin(0,0,2,0,unit = "mm"),
        plot.background = element_blank(),
        panel.spacing = unit(0, "cm"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.grid = element_blank(),
        #strip.text.y.left = element_text(angle=0),
        #strip.background = element_rect(fill="white",linewidth = 0),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(size=10, colour = "black",angle=60, hjust = 1),
        axis.ticks.x= element_line(colour = NA),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        # axis.text.y = element_text(size=8),
        axis.title.y= element_blank()) +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0.5))+
  geom_tile()+ scale_fill_gradient(low="grey95",high="black") +
  labs(title="") +
  xlab("")

library(gridExtra)
lay <- rbind(c(rep(1,33),
               rep(2,1),
               rep(3,5)))

# draw plot
diet.list <- grid.arrange(a,b,c,layout_matrix=lay)
ggsave(filename = "output_files/revised_MAG/Figure3/diet_list.svg",diet.list,device = "svg",width = 8.5, height = 9.5,units = "in")
saveRDS(grid.arrange(a,b,c,layout_matrix=lay),
        file = "output_files/revised_MAG/Figure3/Fig3E.rds")






### POP RANEF ####
# create figure data frame
ggdata <- popranef.df
ggdata$topPop <- factor(ggdata$topPop, c("Asmat_DAI","Basap_BRU","Punan_HUL","Balinese_DPS", "Balinese_PDW","none"))
ggdata <- ggdata[ggdata$pop_var!=0,]
ggdata$label <- ggdata$feature

ggdata <- merge.data.frame(ggdata, full.taxonomy, by="feature")
write.table(ggdata,"output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/Sago+Chicken+Pork_popranef/all_results(Sago+Pork+Chicken)_popvar.tsv", quote = F, row.names = F, sep = "\t")


summary(ggdata$pop_var)
quantile(ggdata$pop_var)
hist(ggdata$pop_var, breaks=30)

ggdata %>% summarise(n=length(pop_var), 
                     average=mean(pop_var),
                     sd=sd(pop_var),
                     se=sd(pop_var)/length(pop_var),
                     min=min(pop_var),
                     Q1=quantile(pop_var)[2],
                     med=median(pop_var),
                     Q3=quantile(pop_var)[4],
                     max=max(pop_var))

ggdata %>% 
  filter(pop_var > 35) %>% 
  group_by(phylum) %>%
  summarise(n=length(pop_var), 
            average=mean(pop_var),
            sd=sd(pop_var),
            se=sd(pop_var)/length(pop_var),
            min=min(pop_var),
            Q1=quantile(pop_var)[2],
            med=median(pop_var),
            Q3=quantile(pop_var)[4],
            max=max(pop_var))



# Filter MAGs 
ggdata <- ggdata[order(ggdata$pop_var, decreasing = T),]
best.pop <- ggdata[ggdata$pop_var>35,] # varex > 35%
# best.pop <- ggdata 

# set phylum colours
phylum.cols <- read.delim("input_files/phylum_colours_bar.tsv", header=T)
tmp <- phylum.cols$cols
names(tmp) <- phylum.cols$phylum
phylum.cols <- tmp

# set population colours
popcols <- readRDS("input_files/pop_mycols.rds")
popcols <- popcols[names(popcols) %in% unique(input.metadata$population)]

table(best.pop$phylum)

# order data
best.pop <- best.pop[order(best.pop$pop_var,decreasing = T),]
best.pop$feature <- factor(best.pop$feature, level=unique(best.pop$feature))
mag.levels <- unique(best.pop$feature)
best.pop$feature <- factor(best.pop$feature, rev(mag.levels))
best.pop <- best.pop[order(rev(mag.levels)),]

phylum.tree.order <- rev(readRDS("output_files/revised_MAG/phylum.tree.order.rds"))
phylum.tree.order <- c(phylum.tree.order,"Methanobacteriota")
best.pop$phylum <- factor(best.pop$phylum, phylum.tree.order)

rownames(best.pop) <- NULL

# add diet associations
diet.assoc <- df[df$qval<0.05, c("feature","value")]
best.pop$diet.assoc <- "population"
for(i in 1:nrow(best.pop)){
  x <- as.character(best.pop$feature)[i]
  y <- x %in% diet.assoc$feature
  z <- diet.assoc$value[diet.assoc$feature==x]
  if(isTRUE(y)){
    best.pop[i,"diet.assoc"] <- y
  }
}

########
# plot #
########
diet.cols <- c("#00A600","#FF69B4","#0000FF","#7A67EE","grey70")
names(diet.cols) <- c("Sago","Pork","Chicken","population","not significant")

#ggplot(best.pop, aes(x=pop_var,y=feature))+ theme_classic() +
#  theme(plot.margin = margin(0,0,0,0,unit = "mm"),
#        panel.spacing = unit(0, "cm"),
#        legend.position = "right",
#        panel.background = element_rect(color = "black",fill=NA),
#        strip.text.y.left = element_text(angle=0),
#        strip.background = element_rect(fill="white",linewidth = 0),
#        text=element_text(size=12), axis.text=element_text(size=12),
#        axis.text.y = element_text(angle=0, hjust = 1, size=8),
#        axis.text.x=element_text(colour = "black", size=8),
#        # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#        axis.title.y= element_blank(),
#        panel.grid.major.x = element_line(colour="grey",linetype = 2)) + 
#  scale_color_manual(values = diet.cols) +
#  scale_x_continuous(breaks = seq(0,100, by=10),limits = c(35,105), expand = c(0,0),labels = paste0(seq(0,100, by=10),"%"))+
#  scale_y_discrete(expand = c(0,0.5))+
#  xlab("Variance Explained") +
#  labs(title="MAGs with Strongest Population Effect (varience explained >35%)")+
#  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
#  geom_point(aes(col=diet.assoc),size=3) +
#  geom_text(aes(label=species), angle=0, nudge_x = 0.5, 
#            size=3, col="grey40", hjust=0, vjust=0.2)

# Total Abundance
mag.levels <- as.character(unique(best.pop$feature))

df.abund <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_comm.tsv",header = T)
df.abund$sampleid <- rownames(df.abund)
df.abund <- melt(df.abund, id.vars = "sampleid",variable.name = "feature",value.name = "CPM")

input.metadata <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_metadata.tsv",header = T)
input.metadata$sampleid <- rownames(input.metadata)
input.metadata <- input.metadata[,c("sampleid","population")]

df.abund <- merge.data.frame(input.metadata,df.abund, by="sampleid")

df.abund <- merge.data.frame(df.abund, full.taxonomy, by="feature")

n <- length(unique(df.abund$sampleid))

df.abund %>% group_by(feature, phylum) %>% summarise(abundance=sum(CPM)/(n*10^6)) -> total_abund

# order data
heat.df.total <- total_abund[total_abund$feature %in% as.character(mag.levels),]
heat.df.total$feature <- factor(heat.df.total$feature, mag.levels)
heat.df.total$phylum <- factor(heat.df.total$phylum, phylum.tree.order)

########
# plot #
########

# ggplot(heat.df.total)+ theme_classic() +
#  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
#  geom_tile(aes(x="%Abundance (log10)", y=feature,fill=log(abundance)))+
#  scale_fill_gradient(low = "grey90", high = "black",)+
#  theme(plot.margin = margin(0,0,0,0,unit = "mm"),
#        panel.spacing = unit(0, "cm"),
#        legend.position = "right",
#        panel.grid = element_blank(),
#        panel.border = element_rect(colour = "black", fill=NA),
#        strip.text.y.left = element_text(angle=0),
#        strip.background = element_rect(fill="white",linewidth = 0),
#        text=element_text(size=12), axis.text=element_text(size=12),
#        axis.ticks.x= element_blank(),
#        axis.text.x=element_text(colour = "black", size=8),
#        # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#        axis.text.y = element_text(size=8),
#        axis.ticks.y= element_blank(),
#        axis.title.y= element_blank()) +
#  labs(title="Total Abundance") +
#  scale_x_discrete(expand = c(0,0))+
#  scale_y_discrete(expand = c(0,0.5))+
#   xlab("")


# add heatmap abundance # normalised average
df.abund <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_comm.tsv",header = T)
df.abund$sampleid <- rownames(df.abund)
df.abund <- melt(df.abund, id.vars = "sampleid",variable.name = "feature",value.name = "CPM")

input.metadata <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_metadata.tsv",header = T)
input.metadata$sampleid <- rownames(input.metadata)
input.metadata <- input.metadata[,c("sampleid","population")]

df.abund <- merge.data.frame(input.metadata,df.abund, by="sampleid")

df.abund <- merge.data.frame(df.abund, full.taxonomy, by="feature")

n <- length(unique(df.abund$sampleid))

df.abund %>% group_by(population,feature, phylum) %>% 
  summarise(popCPM=sum(CPM),
            popMedian=median(CPM),
            popMean=mean(CPM),
            popSD=sd(CPM)) -> pop_abund



df.abund %>% group_by(feature) %>% summarise(totalCPM=sum(CPM)) -> total_abund

heat.df.pop <- merge.data.frame(pop_abund, total_abund, by="feature")

pop_abund %>% group_by(feature) %>% summarise(sum.popMean=sum(popMean)) -> sumpopMean
heat.df.pop <- merge.data.frame(heat.df.pop, sumpopMean, by="feature")

heat.df.pop <- heat.df.pop[order(heat.df.pop$feature),]

heat.df.pop$normalised.abund <- 100*heat.df.pop$popCPM/heat.df.pop$totalCPM
heat.df.pop$normalised.average <- 100*heat.df.pop$popMean/heat.df.pop$sum.popMean

# order data
heat.df.pop <- heat.df.pop[heat.df.pop$feature %in% as.character(mag.levels),]
heat.df.pop$feature <- factor(heat.df.pop$feature,mag.levels)
heat.df.pop$phylum <- factor(heat.df.pop$phylum, phylum.tree.order)
heat.df.pop$population <- factor(heat.df.pop$population,
                                 c("Asmat_DAI","Punan_HUL","Basap_BRU",
                                   "Balinese_PDW","Balinese_DPS"))

########
# plot #
########

# ggplot(heat.df.pop, aes(x=population, y=feature, fill=normalised.average)) + 
#  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
#  theme_minimal() + 
#  theme(plot.margin = margin(0,0,0,0,unit = "mm"),
#        panel.spacing = unit(0, "cm"),
#        legend.position = "right",
#        panel.grid = element_blank(),
#        panel.border = element_rect(colour = "black", fill=NA),
#        strip.text.y.left = element_text(angle=0),
#        strip.background = element_rect(fill="white",linewidth = 0),
#        text=element_text(size=12), axis.text=element_text(size=12),
#        axis.ticks.x= element_blank(),
#        axis.text.x=element_text(colour = "black", size=8),
#        # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#        axis.text.y = element_text(size=8),
#        axis.ticks.y= element_blank(),
#        axis.title.y= element_blank()) +
#  scale_x_discrete(expand = c(0,0))+
#  scale_y_discrete(expand = c(0,0.5))+
#  geom_tile()+ scale_fill_gradient(low="white",high="black") +
#  labs(title="Normalised Average") +
#  xlab("")

# combine graph
a <- ggplot(best.pop, aes(x=pop_var,y=feature,unit = "px"))+ 
  theme_minimal() + 
  theme(plot.margin = margin(0,0,17,0, unit = "mm"),
        plot.background = element_blank(),
        panel.spacing = unit(0, "cm"),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour="grey",linetype = 3, linewidth = .5),
        panel.grid.minor.x = element_line(colour="grey",linetype = 3, linewidth = .5),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill = NA),
        strip.text.y.left = element_text(angle=0),
        strip.background = element_rect(fill="white",linewidth = 0),
        text=element_text(size=12), axis.text=element_text(size=12),
        #axis.text.y = element_text(angle=0, hjust = 1, size=8),
        axis.text.x=element_text(colour = "black", size=8),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y= element_blank())+
  scale_color_manual(values = diet.cols) +
  scale_x_continuous(breaks = seq(0,100, by=10),limits = c(35,105), 
                     expand = c(0.06,0.06),labels = paste0(seq(0,100, by=10),"%"))+
  scale_y_discrete(expand = c(0,0.5))+
  xlab("Variance Explained") +
  labs(title="MAGs with Population Effect (variance explained >35%)") +
  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
  geom_point(aes(col=diet.assoc),size=3)+
  geom_text(aes(label=species), angle=0, nudge_x = 1.2, 
            size=3, col="grey40", hjust=0, vjust=0.5)


b <- ggplot(heat.df.total)+ 
  theme_minimal() + 
  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
  geom_tile(aes(x="Total Abundance", y=feature,fill=log(abundance)))+
  scale_fill_gradient(low = "grey90", high = "black",)+
  theme(plot.margin = margin(0,0,0,0,unit = "mm"),
        plot.background = element_blank(),
        panel.spacing = unit(0, "cm"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        # strip.text.y.left = element_text(angle=0),
        # strip.background = element_rect(fill="white",linewidth = 0),
        strip.text = element_blank(),
        strip.background = element_blank(),
        text=element_text(size=12), axis.text=element_text(size=12),
        axis.ticks.x= element_blank(),
        axis.text.x=element_text(colour = "black", size=8, angle = 60, hjust = 1),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        #axis.text.y = element_text(size=8), axis.ticks.y= element_blank(),
        axis.title.y= element_blank()) +
  labs(title="") +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0.5))+
  xlab("")

c <- ggplot(heat.df.pop, aes(x=population, y=feature, fill=normalised.average)) + 
  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
  theme_minimal() + 
  theme(plot.margin = margin(0,0,3,0,unit = "mm"),
        plot.background = element_blank(),
        panel.spacing = unit(0, "cm"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        # strip.text.y.left = element_text(angle=0),
        # strip.background = element_rect(fill="white",linewidth = 0),
        strip.text = element_blank(),
        strip.background = element_blank(),
        text=element_text(size=12), axis.text=element_text(size=12),
        axis.ticks.x= element_blank(),
        axis.text.x=element_text(colour = "black", size=8,angle = 60, hjust = 1),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        # axis.text.y = element_text(size=8),axis.ticks.y= element_blank(),
        axis.title.y= element_blank()) +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0.5))+
  geom_tile()+ scale_fill_gradient(low="grey95",high="black") +
  labs(title="") +
  xlab("")

lay <- rbind(c(rep(1,33),
               rep(2,1),
               rep(3,4)))
# draw plot
popvar.list <- grid.arrange(a,b,c,layout_matrix=lay)
ggsave(filename = "output_files/revised_MAG/Figure3/popvar_list.svg", popvar.list, device = "svg",width = 8.5, height = 10,units = "in")
saveRDS(popvar.list, "output_files/revised_MAG/Figure3/Fig3F.rds")
