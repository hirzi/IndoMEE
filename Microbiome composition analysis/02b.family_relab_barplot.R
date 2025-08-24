library(dplyr)
library(reshape2)

setwd()

## mypop ##
df <- read.delim("input_files/MAGs/revised/Hybrid_wSPMP/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv", header = T)
mypop <- readRDS("input_files/mypop_wreplicates.rds")
sampleid <- unique(as.character(mypop$sampleid))
pop.ord <- levels(mypop$population)


bio <- readRDS("output_files/clindata_imputed.rds")[,c("sampleid","Age","Gender")]
bio <- merge.data.frame(mypop[,c("sampleid","population","urbanisation")], bio, by="sampleid")
sinpop <- read.delim("input_files/MAGs/PRJEB49168_Illumina_WGS_samples_metadata.txt")[,c("run_accession","age","gender")]
colnames(sinpop) <- c("sampleid","Age","Gender")
sinpop <- merge.data.frame(mypop[,c("sampleid","population","urbanisation")], sinpop, by="sampleid")
mypop <- rbind(bio, sinpop)
mypop <- mypop[order(mypop$sampleid),]
mypop <- mypop[order(mypop$population),]

table(mypop$population)
sum(table(mypop$population))


##################
# Parse Taxonomy #
##################
# source("parse_taxonomy.r")
# full.taxonomy <- parse_taxonomy(df)
# full.taxonomy$species <- deduplicate_species(full.taxonomy$species,type = "numbers")
# write.table(full.taxonomy,"input_files/MAGs/revised/spore_taxonomy.csv",quote = F,sep = ",",row.names = F)

#####################################
#### Taxonomy and Taxonomy Order ####
#####################################
cutoff <- 4

spore.tax <- read.csv("input_files/MAGs/revised/spore_taxonomy.csv")
colnames(spore.tax)[1] <- "feature"

# order taxon my abundance (Using RPK)
df <- read.delim("input_files/MAGs/revised/Hybrid_wSPMP/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_rpkm_revised.tsv")
tmp <- df[,c("original_bin",sampleid)]
colnames(tmp)[1] <- "feature"
tmp <- merge.data.frame(tmp,spore.tax, by="feature")
tmp <- melt(tmp, id.vars = c(colnames(spore.tax)),variable.name = "sampleid",value.name = "RPK")
tmp <- tmp[tmp$RPK>0,]
tmp <- merge.data.frame(tmp, mypop[,c(1,2)], by="sampleid")

### phylum
tmp %>% 
  group_by(phylum,sampleid) %>% summarise(n=length(RPK), total=mean(RPK)) %>%
  group_by(phylum) %>% summarise(n=length(total), average=sum(total)) -> phylum.abund
phylum.abund <- phylum.abund[order(phylum.abund$average, decreasing = T),]
phylum.ord <- unique(phylum.abund$phylum)
phylum.ord <- c(phylum.ord[which(phylum.ord=="Actinomycetota")],phylum.ord[-which(phylum.ord=="Actinomycetota")])
phylum.abund
length(phylum.ord)

### family
tmp %>% 
  group_by(phylum,family,sampleid) %>% summarise(n=length(RPK), total=mean(RPK)) %>%
  group_by(phylum,family) %>% summarise(n=length(total), average=sum(total)) -> family.abund
family.abund <- family.abund[order(family.abund$average, decreasing = T),]
family.abund

### without MAASLIN
family.select.df <- family.abund

### significant taxon
# maaslin.output <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/cont_p20/lifestyle_significant_results_wtax.tsv")
# signif.family <- unique(maaslin.output$family)
# signif.family <- distinct(subset(spore.tax[,c("phylum","family")], family %in% signif.family))
# table(signif.family$phylum)
# family.select.df <- subset(family.abund, family %in% signif.family$family)

## Note that phylum with only 1 family is removed for visual coherence 
selected.family <- c()
for(p in phylum.ord){
  subf <- subset(family.select.df, phylum==p)
  lenf <- nrow(subf)
  
  if(lenf>cutoff){
    top3f <- subf$family[1:cutoff]
    selected.family <- c(selected.family, as.character(top3f))
  } else {
    if(lenf>1) {
      selected.family <- c(selected.family, as.character(subf$family))
    }
  }
}
print(selected.family)

selected.family

family.abund[family.abund$family%in% selected.family,] -> x
y <- data.frame(table(x$phylum))
colnames(y) <- c("phylum","family.count")
y$family <- c()
for(i in y$phylum) {
  z <- paste0(x$family[x$phylum==i], collapse = ",")
  y[y$phylum==i,"family"] <- z
}
# y -> family.list.wMAASLIN
y -> family.list

print(family.list)
 
family.abund0 <- family.abund
family.abund <- subset(family.abund, family %in% selected.family)
family.abund$phylum <- factor(family.abund$phylum,phylum.ord, phylum.abund$phylum)
family.abund <- family.abund[order(family.abund$phylum),]
family.ord <- c(unique(as.character(family.abund$family)),"Others")

reversed.family.abund <-  family.abund[order(family.abund$average),]
reversed.family.abund$phylum <- factor(reversed.family.abund$phylum, levels = rev(phylum.ord))
reversed.family.abund <- reversed.family.abund[order(reversed.family.abund$phylum),]
reversed.family.ord <- c("Others",reversed.family.abund$family)


###############
### barplot ###
###############
library(dplyr)
library(reshape2)

pop.ord <- levels(bio$population)

# prepare data
df <- read.delim("input_files/MAGs/revised/Hybrid_wSPMP/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_rpkm_revised.tsv")
RPK <- df[,sampleid]
rownames(RPK) <- df$original_bin

# remove obese
Ob <- readRDS("input_files/Ob_samples.rds")
RPK <- RPK[!rownames(RPK) %in% Ob,]

RPK <- data.frame(as.matrix(t(RPK)))
RPK$sampleid <- rownames(RPK)
RPK <- melt(RPK,variable.name = "feature",value.name = "RPK")
RPK <- merge.data.frame(RPK, mypop[,c("sampleid","population")], by="sampleid")
RPK <- merge.data.frame(RPK, spore.tax, by="feature")
RPK[,c(5:13)] <- apply(RPK[,c(5:13)], 2, function(x) as.character(x))

RPK %>% group_by(sampleid) %>% summarise(sTotal=sum(RPK)) -> RPK_stotal

RPK %>% filter(family %in% selected.family) %>% select(feature) %>% unique() -> selected.features
selected.features <- as.character(selected.features$feature)

RPK %>% filter(feature %in% selected.features) %>% 
  group_by(sampleid, phylum, family) %>% 
  summarise(sTotal=sum(RPK)) -> selected_family_total

selected_family_total %>% group_by(sampleid) %>% summarise(selecTotal=sum(sTotal)) -> select_total

tmp <- merge.data.frame(RPK_stotal, select_total, by="sampleid")
tmp$Others <- tmp$sTotal-tmp$selecTotal
othersTotal <- tmp[,c("sampleid","Others")]
colnames(othersTotal) <- c("sampleid","sTotal")
othersTotal$family <- "Others"
othersTotal$phylum <- "Others"

selected_family_total <- selected_family_total[order(selected_family_total$sTotal,decreasing = T),]

family.total <- rbind(selected_family_total[,c("sampleid","phylum","family","sTotal")],
                      othersTotal[,c("sampleid","phylum", "family","sTotal")])
colnames(family.total)<-c("sampleid","phylum","family","total")

unique(family.total$family)

family.total$family <- factor(family.total$family, levels = family.ord)

top.family <- family.abund$family[1]
sample.ord <- subset(family.total, family==top.family)
sample.ord <- sample.ord[order(sample.ord$total,decreasing = T),]
sample.ord <- unique(as.character(sample.ord$sampleid))

RA <- recast(family.total[,c("phylum","family","sampleid","total")], family ~ sampleid)
rownames(RA) <- as.character(RA$family)
RA <- RA[,-1]
RA <- apply(RA, 2, function(x) x/sum(x))
RA <- RA[,sample.ord]
RA <- data.frame(as.matrix(t(RA)))
RA$sampleid <- rownames(RA)
RA <- melt(RA,variable.name = "family",value.name = "RA",factorsAsStrings = T)
RA <- merge.data.frame(mypop[,c("sampleid","population")], RA, by="sampleid")
RA$family <- as.character(RA$family)
RA$family <- gsub(pattern = "\\.","-",x = RA$family)
RA$family <- factor(RA$family, levels=family.ord)
RA$sampleid <- factor(RA$sampleid, levels = sample.ord)
RA <- RA[order(RA$sampleid),]
RA <- RA[order(RA$family),]

########################
#### Colour Palette ####
########################
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

# functions
preview_colours <- function(labels, cols){
  barplot(rep(1,length(labels)), col=cols,axes = F, xpd = T, space=0)
  text(data.frame(x=c(1:length(labels))-0.5, y=0.05),labels = labels,srt=90,adj = 0)
}
TaxonColRamp <- function(taxon.data,selected.taxon,level, phylum.cols) {
  a <- grep(as.character(level), colnames(taxon.data))
  a <- as.character(taxon.data[,a])
  b <- taxon.data[a %in% selected.taxon,c("phylum",level)]
  colnames(b)[2] <- level
  b$collate <- paste0(b[,"phylum"],";",b[,level])
  b <- unique(b$collate)
  b <- data.frame(phylum=unlist(lapply(strsplit(b,split = ";"), function(x) x[1])),
                  level=unlist(lapply(strsplit(b,split = ";"), function(x) x[2])))
  colnames(b)[2] <- level
  b <- merge.data.frame(b,phylum.cols, by="phylum")
  
  for(i in unique(b$phylum)){
    myphil <- as.character(b$phylum)
    len <- nrow(b[myphil==i,])
    startcol <- unique(b[myphil==i,"cols"])
    c <- colorRampPalette(c(startcol,"grey90"))(8)
    d <- colorRampPalette(c(startcol,c[6]))(len)
    b[myphil==i,"cols"] <- d
  }
  b <- b[order(b$phylum),]
  
  len <- nrow(b)
  b[len+1,"phylum"] <- "Others"
  b[len+1,level] <- "Others"
  b[len+1,"cols"] <- "grey90"
  b
}
col_df_to_vector <- function(cols,names) {
  tmp <- cols
  names(tmp) <- names
  tmp
}

mycols <- read.delim("input_files/MAGs/MAGs_pop_mycols.tsv", header=T)
preview_colours(labels = mycols$group, cols = mycols$mycols)

phylum.cols <- read.delim("input_files/phylum_colours_bar.tsv", header=T)
phylum.cols[phylum.cols$phylum=="Bacillota","cols"] <- "#DAA520"
phylum.cols[phylum.cols$phylum=="Bacillota_B","cols"] <- "#9ACD32"
phylum.cols[phylum.cols$phylum=="Pseudomonadota","cols"] <- "#473C8B"
phylum.cols[phylum.cols$phylum=="Methanobacteriota","cols"] <- "#CD6090"
phylum.cols[phylum.cols$phylum=="Cyanobacteriota","cols"] <- "#1E90FF"
phylum.cols[phylum.cols$phylum=="Campylobacterota","cols"] <- "#AD9ACA"

preview_colours(labels = phylum.cols$phylum, cols = phylum.cols$cols)

as.character(unique(RA$family)) [!as.character(unique(RA$family)) %in% unique(spore.tax$family)]

taxdata <- distinct(subset(spore.tax[,c("phylum","family")],family%in%selected.family))
taxdata$family <- factor(taxdata$family, family.ord)
taxdata <- taxdata[order(taxdata$family),]
family.cols <- TaxonColRamp(taxon.data = taxdata,
                            selected.taxon = selected.family,
                            level = "family",
                            phylum.cols = phylum.cols)
family.cols$family <- factor(family.cols$family, levels = family.ord)
family.cols <- family.cols[order(family.cols$family),]
row.names(family.cols) <- NULL
preview_colours(labels = family.cols$family, cols = family.cols$cols)

legend2D <- ggplot(family.cols) + 
                    theme_void(base_size = 8) +
                    geom_bar(aes(x=1, y=1, fill=family),position = "stack", stat = "identity")+
                    geom_text(aes(x=1, y=-0.5+1:nrow(family.cols), label=rev(family)))+
                    scale_fill_manual(values = family.cols$cols) + 
                    theme(plot.margin = margin(0,0,0,0,"mm"), legend.position = "none")

legend2D
saveRDS(legend2D,
        file = "output_files/revised_MAG/Figure2/Fig2D.legend.rds")
saveRDS(family.cols, file = "input_files/family_cols.rds")

# convert to vector
mycols <- col_df_to_vector(cols = mycols$mycols, names = mycols$group)
phylum.cols <- col_df_to_vector(cols = phylum.cols$cols, names = phylum.cols$phylum)
family.cols <- col_df_to_vector(cols = family.cols$cols, names = family.cols$family)


### Barplot Data
ggdata <- RA
cluster_df <- read.delim("output_files/revised_MAG/ward_clusters.tsv")
ggdata <- merge.data.frame(ggdata, cluster_df[,c("sampleid","clust_k3")], by="sampleid")
ggdata$cluster <- ifelse(ggdata$clust_k3==1, "Cluster 1",
                         ifelse(ggdata$clust_k3==2, "Cluster 3", "Cluster 2"))

# ward cluster sample sort
# library(ape)
# wardtree <- read.tree("output_files/MAGs/SIN_tree.trefile")
# node.ord <- wardtree$tip.label
# ggdata$sampleid <- factor(ggdata$sampleid, node.ord) # classic sort

check.completeness <- as.character(unique(ggdata$family))
table(check.completeness %in% as.character(family.ord))
table(check.completeness %in% names(family.cols))
check.completeness[!check.completeness %in% family.ord]
check.completeness[!check.completeness %in% names(family.cols)]

ggdata$family <- factor(ggdata$family, reversed.family.ord)
ggdata$sampleid <- factor(ggdata$sampleid, rev(sample.ord))
ggdata$population <- factor(ggdata$population, pop.ord)
rownames(ggdata) <- NULL

# panel A
par(mar=c(3,3,3,3))
a <- ggplot(ggdata, 
       aes(x=sampleid, y=RA, fill=family, group=family))+
  theme_minimal()+
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        axis.text= element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey40", fill=NA),
        panel.spacing = unit(0.1, "cm", data = NULL),
        strip.text=element_text(angle=90)) +
  facet_grid(~population,scale = "free",space = "free")+
  geom_bar(position = "stack", stat="identity", width=1) +
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(expand=c(0,0)) +
  scale_fill_manual(values=family.cols) +
  guides(fill = guide_legend(reverse=TRUE,nrow = cutoff)) 

#  ggtitle(paste("Top",cutoff,"Family in each Phylum housing a significant MAG"))

a

## by cluster
clustercols <- c("grey80","grey50","grey30")
names(clustercols) <- c("Cluster 1", "Cluster 2", "Cluster 3")

b <-ggplot(ggdata, 
           aes(x=sampleid, y=1))+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey40", fill=NA),
        panel.spacing = unit(0.1, "cm", data = NULL),
        strip.text = element_text(angle=90)) +
  facet_grid(~population,scale = "free",space = "free")+
  geom_tile(aes(fill=cluster)) +
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(expand=c(0,0)) +
  scale_fill_manual(values=clustercols)

b

lay <- rbind(2,1,1,1)
grid.arrange(a+theme(strip.text = element_blank()),b,layout_matrix=lay)

ggsave(grid.arrange(a+theme(strip.text = element_blank()),b,layout_matrix=lay),
       filename = "output_files/revised_MAG/Figure2/family_bar_cluster_annotated.svg",
       height = 16, width = 7,dpi = 300)

saveRDS(a, "output_files/revised_MAG/Figure2/Fig2D.2.rds")
saveRDS(b, "output_files/revised_MAG/Figure2/Fig2D.1.rds")


a1 <- a+theme(strip.text = element_blank(),
             legend.position = "none",
             plot.margin = margin(0,0,0,0.5,unit = "cm"))

b <- b + theme(plot.margin = margin(0,0,0,0.5,unit = "cm"))


lay <- rbind(2,2,1,1,1,1)
c <- grid.arrange(a1,
                  b,
                  layout_matrix=lay)
saveRDS(c,"output_files/revised_MAG/Figure2/Fig2D.collate.rds" )






##################
# compare family #
##################
#clust_df_noRPN <- cluster_df[cluster_df$population!="Punan_RPN",]
clust_df_noRPN <- cluster_df
chisq.test(table(clust_df_noRPN$urbanisation, clust_df_noRPN$clust_k3))

family.total <- RPK

family.total %>% 
  group_by(phylum, family, sampleid, population) %>% 
  summarise(total=sum(RPK)) -> family.total

family.levels <- unique(family.total$family,"Others")
family.total$family <- factor(family.total$family, levels = family.levels)

phylum.levels <- unique(family.total$phylum, "Others")
family.total$phylum <- factor(family.total$phylum,phylum.levels)

# sample.ord <- family.total[family.total$family==selected.family[1],]
sample.ord <- family.total[,]
sample.ord <- sample.ord[order(sample.ord$total,decreasing = F),]
sample.ord$population <- factor(sample.ord$population,pop.ord)
sample.ord <- sample.ord[order(sample.ord$population),]
sample.ord <- unique(sample.ord$sampleid)

family.total <- recast(family.total[,c("family","sampleid","total")], family ~ sampleid)
rownames(family.total) <- family.total$family
family.total <- family.total[,-1]
family.total <- apply(family.total,2, function(x) x/sum(x))
family.total <- family.total[,sample.ord]

RA <- apply(family.total,2, function(x) x/sum(x))
RA <- data.frame(family=rownames(RA), RA[,sample.ord])

RA <- melt(RA, id.vars="family", variable.name = "sampleid", value.name = "RA")

kruskal.df <- merge.data.frame(RA, cluster_df[,c("sampleid","population","clust_k3","urbanisation")], by="sampleid")
kruskal.df$cluster <- ifelse(kruskal.df$clust_k3==1, "Cluster 1",
                             ifelse(kruskal.df$clust_k3==3, "Cluster 3", "Cluster 2"))

kruskal.out <- list()
wilcox.out <- list()
effect.out <- list()
for(i in unique(kruskal.df$family)){
  tmp.df <- subset(kruskal.df, family==i)
  krus.out <- kruskal.test(RA~cluster, data=tmp.df)
  pwc.out <- pairwise.wilcox.test(tmp.df$RA, g=tmp.df$cluster, paired = F,p.adjust.method = "BH")
  eff.out <- rstatix::kruskal_effsize(RA~cluster, data=tmp.df,ci = T,nboot = 100)
  out <- list(Kruskal=list(krus.out), Wilcox=list(pwc.out), Effect=eff.out$effsize[1])
  n <- names(kruskal.out)
  n <- c(n,i)
  len <- length(kruskal.out)
  kruskal.out[len+1] <- list(krus.out)
  names(kruskal.out) <- n
  wilcox.out[len+1] <- list(pwc.out)
  names(wilcox.out) <- n
  effect.out[len+1] <- list(eff.out)
  names(wilcox.out) <- n
  }

krus.signif <- data.frame(family=names(kruskal.out),
                          stat=unlist(lapply(kruskal.out, function(x) x$statistic)),
                          pval=unlist(lapply(kruskal.out, function(x) x$p.value)), 
                          effsize=unlist(lapply(effect.out, function(x) x$effsize)),
                          conf.low=unlist(lapply(effect.out, function(x) x$conf.low)),
                          conf.high=unlist(lapply(effect.out, function(x) x$conf.high)),
                          row.names = NULL)
krus.signif$qval <- p.adjust(krus.signif$pval)
write.table(krus.signif, file = "output_files/revised_MAG/SPMP_MAGS_Kruskal_Family_vs_Cluster.tsv", row.names = F, quote = F, sep = "\t")
table(krus.signif$qval<0.05)
krus.signif <- krus.signif[krus.signif$qval < 0.05, ]
krus.signif <- krus.signif[order(krus.signif$stat, decreasing = T),]
View(krus.signif)

signif.fam <- krus.signif$family

wilcox.signif <- data.frame()
for(i in unique(RA$family)){
  out <- wilcox.out[[i]]$p.value
  wilcox.tmp <- data.frame(family=i,
                           melt(out))
  wilcox.signif <- rbind(wilcox.signif, wilcox.tmp)
}
wilcox.signif <- na.omit(wilcox.signif)
table(wilcox.signif$value<0.05)
wilcox.signif <- wilcox.signif[wilcox.signif$value < 0.05,]
wilcox.signif <- wilcox.signif[order(wilcox.signif$Var1,decreasing = T),]
write.table(krus.signif, file = "output_files/revised_MAG/SPMP_MAGS_Wilcox_Family_vs_Cluster_Pair.tsv", row.names = F, quote = F, sep = "\t")

kruskal.df %>% 
  group_by(population) %>% 
  summarise(maxiC0=max(RA), mediC0=median(RA), aver=mean(RA)) -> maxi_medi_bypop
maxi_medi_bypop

kruskal.df$vsCluster3 <- ifelse(kruskal.df$cluster=="Cluster 3", "Cluster 3", "Cluster 1&2")

RA_subs <- melt(family.total,varnames = c("family","sampleid"))
RA_subs$value <- RA_subs$value + 10^-5
RA_subs <- merge.data.frame(RA_subs, distinct(kruskal.df[,c("sampleid","vsCluster3")]))
RA_subs %>% group_by(family, vsCluster3) %>% summarise(meanRPK=mean(value, na.rm = T),
                                                       medianRPK=median(value, na.rm = T)) -> RA_subs
mean_RAdiff <- RA_subs[,c("family","vsCluster3","meanRPK")]
medi_RAdiff <- RA_subs[,c("family","vsCluster3","medianRPK")]
mean_RAdiff <- recast(formula=family ~ vsCluster3, mean_RAdiff)
medi_RAdiff <- recast(formula=family ~ vsCluster3, medi_RAdiff)
mean_RAdiff$diff <- mean_RAdiff$`Cluster 3` - mean_RAdiff$`Cluster 1&2`
medi_RAdiff$diff <- medi_RAdiff$`Cluster 3` - medi_RAdiff$`Cluster 1&2`
mean_RAdiff$fold.increase <- c()
medi_RAdiff$fold.increase <- c()
for(i in 1:nrow(mean_RAdiff)){
  mean_RAdiff[i,"fold.increase"] <- max(c(mean_RAdiff[i,2],mean_RAdiff[i,3]))/min(c(mean_RAdiff[i,2],mean_RAdiff[i,3]))
  medi_RAdiff[i,"fold.increase"] <- max(c(medi_RAdiff[i,2],medi_RAdiff[i,3]))/min(c(medi_RAdiff[i,2],medi_RAdiff[i,3]))}
mean_RAdiff$fold.increase <- sign(mean_RAdiff$diff)*round(mean_RAdiff$fold.increase,1)
medi_RAdiff$fold.increase <- sign(medi_RAdiff$diff)*round(medi_RAdiff$fold.increase,1)
colnames(mean_RAdiff)[5] <- "mean.fold.increase"
colnames(medi_RAdiff)[5] <- "median.fold.increase"
merge.data.frame(
  merge.data.frame(
    merge.data.frame(mean_RAdiff,medi_RAdiff[,c(1,5)], by="family"),
                 krus.signif, by="family"),
                 taxdata, by="family") -> maxi_medi

View(subset(maxi_medi, family %in% selected.family))
View(subset(maxi_medi, qval < 0.05))

write.table(maxi_medi,"output_files/revised_MAG/family.kruskal.result.tsv", sep = "\t",quote = F,row.names = F)
write.table(subset(maxi_medi, family %in% selected.family),"output_files/revised_MAG/family.kruskal.result.selected.tsv", sep = "\t",quote = F,row.names = F)
