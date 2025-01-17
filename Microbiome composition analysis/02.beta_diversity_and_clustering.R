# load libraries
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(reshape2)
library(dplyr)
library(vegan)
source("pretty_plot_function.R")

# custom functions
convert_to_CPM <- function(df){
  df <- df*10^6
  df}
prevalence_filter <- function(df, n) {
  x <- df # matrix of abundance with sampleid as columns
  y <- x
  y[y>0] <- 1
  sel <- rowSums(y) 
  print(paste(sum(sel>2), "MAGs are >", n, "prevalence."))
  y <- y[sel>2,]
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
axis_labels <- function(cap){
  eig <- cap$CA$eig
  var <- eig/sum(eig)
  axis_labs <- paste("PCo",1:length(names(var)), " (", round(var*100,1),"%)", sep="")
  axis_labs
}

# load data
df <- read.delim("input_files/MAGs/revised/Hybrid_wSPMP/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv", header = T)
mypop <- readRDS("input_files/mypop_wreplicates.rds")
sampleid <- unique(as.character(mypop$sampleid))
pop.ord <- levels(mypop$population)

# Principal Coordinate Analysis
# With SPMP
# input comm
df <- read.delim("input_files/MAGs/revised/Hybrid_wSPMP/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv", header = T)
sampleid <- readRDS("output_files/sampleid_n116.rds")
sporeids <- readRDS("output_files/sampleid_spore.rds")
reads <- df[,colnames(df) %in% c(sampleid,sporeids)]
rownames(reads) <- df$Genome
reads <- convert_to_CPM(reads)
reads <- prevalence_filter(reads, n=2)
reads <- as.matrix(t(reads))

# input metadata
env.df <- readRDS("input_files/mypop_wSPMP.rds")
env.df <- env.df[env.df$sampleid %in% c(sampleid,sporeids),]
rownames(env.df) <- env.df$sampleid
env.df <- env.df[rownames(reads),]

# Create Ordination Object
# Using Aitchison env.df.imputed# Using Aitchison Distance 
# Aitchison dist transform table with clr
# calculate distance
comm.df <- reads
comm.dist <- vegdist(comm.df+1, method="aitchison")
cap <- capscale(comm.dist~1, data=env.df, comm = comm.df)

# set colours
mycols <- read.delim("input_files/MAGs/MAGs_pop_mycols.tsv", header=T)
mycols <- mycols[mycols$group %in% unique(env.df$population),]
col.ord <- mycols$group
mycols <- mycols$mycols
names(mycols) <- col.ord
mycols <- mycols[pop.ord]

# set shape
myshape <- c(21,22,24)
names(myshape) <- c("remote","rural","urban")

# Create Ordination Object
# Using Aitchison env.df.imputed# Using Aitchison Distance 
# Aitchison dist transform table with clr
# calculate distance
comm.df <- reads
comm.dist <- vegdist(comm.df+1, method="aitchison")
cap <- capscale(comm.dist~1, data=env.df, comm = comm.df)

figpath="output_files/revised_MAG/PCoA/"
figname="revised_wSPMP_PCo1_PCo2"
ggsave(plot = pretty_PCA(cap, axes=c(1,2), metadata = env.df, 
                         colour.order = pop.ord, colour.by = "population", colour.scheme = mycols,
                         shape.by = "urbanisation", shape.order = names(myshape),shape.scheme = myshape),
       filename = paste0(figpath,figname,".png"), 
       device = "png", scale=4, dpi = 600, width = 5, height = 5, units = "cm")


figname="revised_wSPMP_PCo1_PCo3"
ggsave(plot = pretty_PCA(cap, axes=c(1,3), metadata = env.df, 
                         colour.order = pop.ord, colour.by = "population", colour.scheme = mycols,
                         shape.by = "urbanisation", shape.order = names(myshape),shape.scheme = myshape),
       filename = paste0(figpath,figname,".png"), 
       device = "png", scale=4, dpi = 600, width = 5, height = 5, units = "cm")

figname="revised_wSPMP_PCo2_PCo3"
ggsave(plot = pretty_PCA(cap, axes=c(2,3), metadata = env.df, 
                         colour.order = pop.ord, colour.by = "population", colour.scheme = mycols,
                         shape.by = "urbanisation", shape.order = names(myshape),shape.scheme = myshape),
       filename = paste0(figpath,figname,".png"), 
       device = "png", scale=4, dpi = 600, width = 5, height = 5, units = "cm")


# With SPMP on UHGG
# input comm
df <- read.delim("input_files/MAGs/revised/UHGG_wSPMP/bwa_counts_total_filtered_wMetadata_UHGG_strongFilter_relAbund_revised.tsv", header = T)
sampleid <- readRDS("output_files/sampleid_n116.rds")
sporeids <- readRDS("output_files/sampleid_spore.rds")
reads <- df[,colnames(df) %in% c(sampleid,sporeids)]
rownames(reads) <- df$
reads <- convert_to_CPM(reads)
reads <- as.matrix(t(reads))
prevalence_filter(reads, n=2)

# input metadata
env.df <- readRDS("input_files/mypop_wSPMP.rds")
env.df <- env.df[env.df$sampleid %in% c(sampleid,sporeids),]
rownames(env.df) <- env.df$sampleid
env.df <- env.df[rownames(reads),]

# Create Ordination Object
# Using Aitchison env.df.imputed# Using Aitchison Distance 
# Aitchison dist transform table with clr
# calculate distance
comm.df <- reads
comm.dist <- vegdist(comm.df+1, method="aitchison")
cap <- capscale(comm.dist~1, data=env.df, comm = comm.df)


figpath="output_files/revised_MAG/PCoA/"
figname="revised_wSPMP_onUHGG_PCo1_PCo2"
ggsave(plot = pretty_PCA(cap, axes=c(1,2), metadata = env.df, 
                         colour.order = pop.ord, colour.by = "population", colour.scheme = mycols,
                         shape.by = "urbanisation", shape.order = names(myshape),shape.scheme = myshape),
       filename = paste0(figpath,figname,".png"), 
       device = "png", scale=4, dpi = 600, width = 5, height = 5, units = "cm")

figname="revised_wSPMP_onUHGG_PCo1_PCo3"
ggsave(plot = pretty_PCA(cap, axes=c(1,3), metadata = env.df, 
                         colour.order = pop.ord, colour.by = "population", colour.scheme = mycols,
                         shape.by = "urbanisation", shape.order = names(myshape),shape.scheme = myshape),
       filename = paste0(figpath,figname,".png"), 
       device = "png", scale=4, dpi = 600, width = 5, height = 5, units = "cm")

figname="revised_wSPMP_onUHGG_PCo2_PCo3"
ggsave(plot = pretty_PCA(cap, axes=c(2,3), metadata = env.df, 
                         colour.order = pop.ord, colour.by = "population", colour.scheme = mycols,
                         shape.by = "urbanisation", shape.order = names(myshape),shape.scheme = myshape),
       filename = paste0(figpath,figname,".png"), 
       device = "png", scale=4, dpi = 600, width = 5, height = 5, units = "cm")


# IndoMEE only
# input comm
df <- read.delim("input_files/MAGs/revised/Hybrid/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv", header = T)
reads <- df[,colnames(df)%in%sampleid]
rownames(reads) <- df$Genome
reads <- convert_to_CPM(reads)
reads <- prevalence_filter(reads,2)
reads <- as.matrix(t(reads))

# input metadata
env.df <- read.delim("input_files/Punan_160_full_ClinLs_metadata_imputed_wRep.tsv")
env.df <- env.df[env.df$sampleid %in% rownames(reads),]
rownames(env.df) <- env.df$sampleid
env.df <- env.df[rownames(reads),]

# Create Ordination Object
# Using Aitchison env.df.imputed# Using Aitchison Distance 
# Aitchison dist transform table with clr
# calculate distance
comm.df <- reads
comm.dist <- vegdist(comm.df+1, method="aitchison")
cap <- capscale(comm.dist~1, data=env.df, comm = comm.df)

# set colours
mycols <- read.delim("input_files/MAGs/MAGs_pop_mycols.tsv", header=T)
col.ord <- mycols$group
mycols <- mycols$mycols
names(mycols) <- col.ord
mycols <- mycols[pop.ord]

figpath="output_files/revised_MAG/PCoA/"
figname="revised_IndoMEE_PCo1_PCo2"
ggsave(plot = pretty_PCA(cap, axes=c(1,2), metadata = env.df, 
                         colour.order = pop.ord, colour.by = "population", colour.scheme = mycols,
                         shape.by = "urbanisation", shape.order = names(myshape),shape.scheme = myshape),
       filename = paste0(figpath,figname,".png"), 
       device = "png", scale=4, dpi = 600, width = 5, height = 5, units = "cm")

figname="revised_IndoMEE_PCo1_PCo3"
ggsave(plot = pretty_PCA(cap, axes=c(1,3), metadata = env.df, 
                         colour.order = pop.ord, colour.by = "population", colour.scheme = mycols,
                         shape.by = "urbanisation", shape.order = names(myshape),shape.scheme = myshape),
       filename = paste0(figpath,figname,".png"), 
       device = "png", scale=4, dpi = 600, width = 5, height = 5, units = "cm")

figname="revised_IndoMEE_PCo2_PCo3"
ggsave(plot = pretty_PCA(cap, axes=c(2,3), metadata = env.df, 
                         colour.order = pop.ord, colour.by = "population", colour.scheme = mycols,
                         shape.by = "urbanisation", shape.order = names(myshape),shape.scheme = myshape),
       filename = paste0(figpath,figname,".png"), 
       device = "png", scale=4, dpi = 600, width = 5, height = 5, units = "cm")


#################
### Figure 2B ###
#################
df <- read.delim("input_files/MAGs/revised/Hybrid_wSPMP/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv", header = T)
mypop <- readRDS("input_files/mypop_wreplicates.rds")
sampleid <- unique(as.character(mypop$sampleid))
pop.ord <- levels(mypop$population)

# Principal Coordinate Analysis
# With SPMP
# input comm
df <- read.delim("input_files/MAGs/revised/Hybrid_wSPMP/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv", header = T)
sampleid <- readRDS("output_files/sampleid_n116.rds")
sporeids <- readRDS("output_files/sampleid_spore.rds")
reads <- df[,colnames(df) %in% c(sampleid,sporeids)]
rownames(reads) <- df$original_bin
reads <- convert_to_CPM(reads)
reads <- prevalence_filter(reads,n = 2)
reads <- as.matrix(t(reads))

# input metadata
env.df <- readRDS("input_files/mypop_wSPMP.rds")
env.df <- env.df[env.df$sampleid %in% c(sampleid,sporeids),]
rownames(env.df) <- env.df$sampleid
env.df <- env.df[rownames(reads),]

# Create Ordination Object
# Using Aitchison env.df.imputed# Using Aitchison Distance 
# Aitchison dist transform table with clr
# calculate distance
comm.df <- reads
comm.dist <- vegdist(comm.df+1, method="aitchison")
cap <- capscale(comm.dist~1, data=env.df, comm = comm.df)

# set colours
mycols <- read.delim("input_files/MAGs/MAGs_pop_mycols.tsv", header=T)
mycols <- mycols[mycols$group %in% unique(env.df$population),]
col.ord <- mycols$group
mycols <- mycols$mycols
names(mycols) <- col.ord
mycols <- mycols[pop.ord]

p <- pretty_PCA(cap, axes=c(1,2), metadata = env.df, 
                colour.order = pop.ord, colour.by = "population", colour.scheme = mycols,
                shape.by = "urbanisation", shape.order = names(myshape),shape.scheme = myshape)

saveRDS(p, "output_files/revised_MAG/Figure2/Fig2B.rds")
p <- readRDS("output_files/revised_MAG/Figure2/Fig2B.rds")
plot(p)

ggsave(pretty_PCA(cap, axes=c(1,2), metadata = env.df, 
                  colour.order = pop.ord, colour.by = "population", colour.scheme = mycols,
                  shape.by = "urbanisation", shape.order = names(myshape),shape.scheme = myshape),
       filename = "output_files/revised_MAG/Figure2/Fig2B.png", 
       device = "png", scale=4, dpi = 600, width = 5, height = 5, units = "cm")

ggsave(pretty_PCA(cap, axes=c(1,2), metadata = env.df, 
                  colour.order = pop.ord, colour.by = "population", colour.scheme = mycols,
                  shape.by = "urbanisation", shape.order = names(myshape),shape.scheme = myshape),
       filename = "output_files/revised_MAG/Figure2/Fig2B.svg", 
       device = "svg", scale=5, dpi = 600, width = 5, height = 5, units = "cm")


##################
# Save PCoA data #
##################
# SPMP
df <- read.delim("input_files/MAGs/revised/Hybrid_wSPMP/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv", header = T)
sampleid <- readRDS("output_files/sampleid_n116.rds")
sporeids <- readRDS("output_files/sampleid_spore.rds")
reads <- df[,colnames(df) %in% c(sampleid,sporeids)]
rownames(reads) <- df$Genome
reads <- convert_to_CPM(reads)
reads <- prevalence_filter(reads, n=2)
reads <- as.matrix(t(reads))
env.df <- readRDS("input_files/mypop_wSPMP.rds")
env.df <- env.df[env.df$sampleid %in% c(sampleid,sporeids),]
rownames(env.df) <- env.df$sampleid
env.df <- env.df[rownames(reads),]

saveRDS(reads, "output_files/revised_MAG/pmv.comm.MAGs.indospmp.rds")
samples <- rownames(reads)
saveRDS(samples, "output_files/revised_MAG/pmv.sample.indospmp.rds")
saveRDS(env.df, "output_files/revised_MAG/pmv.metadata.indospmp.rds")

# IndoMEE
df <- read.delim("input_files/MAGs/revised/Hybrid/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv", header = T)
sampleid <- readRDS("output_files/sampleid_n116.rds")
reads <- df[,colnames(df) %in% sampleid]
rownames(reads) <- df$Genome
reads <- convert_to_CPM(reads)
reads <- prevalence_filter(reads, n=2)
reads <- as.matrix(t(reads))
env.df <- readRDS("input_files/mypop_wSPMP.rds")
env.df <- env.df[env.df$sampleid %in% c(sampleid,sporeids),]
rownames(env.df) <- env.df$sampleid
env.df <- env.df[rownames(reads),]

saveRDS(reads, "output_files/revised_MAG/pmv.comm.MAGs.indomee.rds")
samples <- rownames(reads)
saveRDS(samples, "output_files/revised_MAG/pmv.sample.indomee.rds")
saveRDS(env.df, "output_files/revised_MAG/pmv.metadata.indomee.rds")

##########################
## Population Centroids ##
##########################
find_centroids <- function(ord, groups, display="sites", axes=c(1,2)){
  centroids <- c()
  for(g in unique(groups)){
    out <- seq(along = groups)
    gr <- out[groups == g]
    sites <- scores(ord, display, choices = axes)
    sites <- sites[gr,]
    centr=t(apply(sites, 2, function(x) weighted.mean(x)))
    rownames(centr) <- g
    centroids <- rbind(centroids, centr)
  }
  centroids
}

# IndoMEE
comm.df <- readRDS("output_files/revised_MAG/pmv.comm.MAGs.indomee.rds")
comm.dist <- vegdist(comm.df+1, method="aitchison")
env.df <- readRDS("output_files/revised_MAG/pmv.metadata.indomee.rds")
cap <- capscale(comm.dist~1, data=env.df, comm = comm.df)

centroids_Indomee <- find_centroids(cap, env.df$population)

# SPMP
comm.df <- readRDS("output_files/revised_MAG/pmv.comm.MAGs.indospmp.rds")
comm.dist <- vegdist(comm.df+1, method="aitchison")
env.df <- readRDS("output_files/revised_MAG/pmv.metadata.indospmp.rds")
cap <- capscale(comm.dist~1, data=env.df, comm = comm.df)
centroids_IndoSPMP <- find_centroids(cap, env.df$population)


###############
## PERMANOVA ##
###############
# IndoMEE
comm.df <- readRDS("output_files/revised_MAG/pmv.comm.MAGs.indomee.rds")
comm.dist <- vegdist(comm.df+1, method="aitchison")
env.df <- readRDS("output_files/revised_MAG/pmv.metadata.indomee.rds")
env.df$cLifestyle <- ifelse(env.df$population %in% c("Asmat_DAI","Punan_PBS","Punan_HUL","Punan_SUL"), "remote",
                            ifelse(env.df$population %in% c("Punan_RPN","Basap_BRU","Balinese_PDW"), "rural",
                                   ifelse(env.df$population %in% c("Balinese_DPS","Chinese_SIN","Malay_SIN","Indian_SIN"),"urban", NA)))
env.df$lifestyle <- ifelse(env.df$population %in% c("Asmat_DAI","Punan_PBS","Punan_HUL","Punan_SUL"), 0,
                           ifelse(env.df$population %in% c("Punan_RPN","Basap_BRU","Balinese_PDW"), 0.5,
                                  ifelse(env.df$population %in% c("Balinese_DPS","Chinese_SIN","Malay_SIN","Indian_SIN"),1, NA)))
sink("output_files/revised_MAG/PERMANOVA_IndoMEE.txt")
adonis2(comm.dist ~ cLifestyle + population, env.df, permutations = 9999)
adonis2(comm.dist ~ lifestyle + population, env.df, permutations = 9999)
sink()

# IndoMEE
comm.df <- readRDS("output_files/revised_MAG/pmv.comm.MAGs.indospmp.rds")
comm.dist <- vegdist(comm.df+1, method="aitchison")
env.df <- readRDS("output_files/revised_MAG/pmv.metadata.indospmp.rds")
env.df$cLifestyle <- ifelse(env.df$population %in% c("Asmat_DAI","Punan_PBS","Punan_HUL","Punan_SUL"), "remote",
                            ifelse(env.df$population %in% c("Punan_RPN","Basap_BRU","Balinese_PDW"), "rural",
                                   ifelse(env.df$population %in% c("Balinese_DPS","Chinese_SIN","Malay_SIN","Indian_SIN"),"urban", NA)))
env.df$lifestyle <- ifelse(env.df$population %in% c("Asmat_DAI","Punan_PBS","Punan_HUL","Punan_SUL"), 0,
                           ifelse(env.df$population %in% c("Punan_RPN","Basap_BRU","Balinese_PDW"), 0.5,
                                  ifelse(env.df$population %in% c("Balinese_DPS","Chinese_SIN","Malay_SIN","Indian_SIN"),1, NA)))
sink("output_files/revised_MAG/PERMANOVA_IndoSPMP.txt")
adonis2(comm.dist ~ cLifestyle + population, env.df, permutations = 9999)
adonis2(comm.dist ~ lifestyle + population, env.df, permutations = 9999)
sink()

##########################
### Sample Clustering ####
##########################
library(ape)

comm.df <- readRDS("output_files/revised_MAG/pmv.comm.MAGs.indospmp.rds")
comm.dist <- vegdist(comm.df+1, method="aitchison")
mypop <- readRDS("output_files/revised_MAG/pmv.metadata.indospmp.rds")
mypop$urbanisation <- factor(mypop$urbanisation, c("remote","rural","urban"))

library(ggdendro)
mycluster <- hclust(comm.dist,method = "ward.D2")
hcdata <- dendro_data(mycluster, type="rectangle")
tipsdata <- label(hcdata)
tipsdata$sampleid <- tipsdata$label
tipsdata <- merge.data.frame(tipsdata, mypop[,c("sampleid","population","urbanisation")])

ward.clust.plot <- ggplot() + 
  theme_void() +
  geom_segment(data=segment(hcdata), aes(x, y, xend=xend, yend=yend))+
  geom_point(data = tipsdata, 
            aes(x = x, y = y,col=population, fill=population), size=1.4, pch=22) +
  scale_colour_manual(values = mycols)+
  scale_fill_manual(values=mycols)
ward.clust.plot
ggsave(plot = ward.clust.plot,"output_files/revised_MAG/Figure2/Fig2C_phylo_SIN.png",device = "png",dpi = 600,width = 15, height = 2,units = "cm",bg = NULL, scale = 2.5)
saveRDS(ward.clust.plot, "output_files/revised_MAG/Figure2/Fig2C_phylo_SIN.rds")
saveRDS(mycluster, "output_files/revised_MAG/ward_clusters.rds")
myclust <- as.phylo(mycluster)
write.tree(myclust,file = "output_files/revised_MAG/SIN_tree.trefile",digits = 5)


# Cluster stats
clust_k4 <- cutree(mycluster,k=4)
clust_k3 <- cutree(mycluster,k=3)
clust_k3[clust_k3==3] <- 4
clust_k3[clust_k3==2] <- 3
clust_k3[clust_k3==4] <- 2
clust_k2 <- cutree(mycluster,k=2)

cluster.df <- data.frame(mypop[mypop$sampleid%in%names(clust_k4),],
                         cbind(clust_k4,
                               clust_k3,
                               clust_k2))
chisq.test(table(cluster.df$urbanisation, cluster.df$clust_k2))
chisq.test(table(cluster.df$urbanisation, cluster.df$clust_k3))
chisq.test(table(cluster.df$urbanisation, cluster.df$clust_k4))

p <- ggplot(cluster.df, aes(x=urbanisation, y=clust_k3, fill=population)) + theme_classic2(base_size = 12) +
  theme(panel.border = element_rect(fill = NA))+
  geom_vline(xintercept = c(1.5,2.5))+
  geom_hline(yintercept = c(1.5,2.5))+
  geom_point(size=4, stat = "identity", pch=21, position = position_jitter(width = 0.45)) +
  scale_fill_manual(values = mycols)
p

write.table(cluster.df, "output_files/revised_MAG/ward_clusters.tsv", row.names = F, quote = F, sep = "\t")
ggsave(p, filename = "output_files/revised_MAG/Figure2/cluster_vs_urbanisation.svg",device = "png", height = 5, width = 5, units = "cm",dpi = 600)

q <- ggplot(cluster.df, aes(x=population, y=clust_k3, fill=population)) + theme_classic2(base_size = 12) +
  theme(panel.border = element_rect(fill = NA))+
  geom_vline(xintercept = seq(1.5,10.5,1))+
  geom_hline(yintercept = c(1.5,2.5))+
  geom_point(size=4, stat = "identity", pch=21, position = position_jitter(width = 0.45)) +
  scale_fill_manual(values = mycols)
q


#k3
tab <- ftable(cluster.df$population,cluster.df$clust_k3)
tab
chisq.test(tab,correct = T,rescale.p = T) #,B = 99999,simulate.p.value = T)
par(mar=c(3,3,3,7), xpd=TRUE)
barplot(as.matrix(tab),ylim = c(0,150),
        col = mycols[sort(unique(cluster.df$population))], 
        main="Count of Samples by Cluster")
legend('right', inset=c(-0.4, 0),pch=15,cex=0.7,
       legend=sort(unique(cluster.df$population)), 
       col=mycols[sort(unique(cluster.df$population))])
tab <- data.frame(tab)
colnames(tab) <- c("Pop","Cluster","Freq")
tab$Percent <- round(100*tab$Freq/sum(tab$Freq),1)
for(i in 1:length(unique(tab$Cluster))){
  mypie <- tab[tab$Cluster==i,]
  mypie$cols <- mycols[as.character(mypie$Pop)]
  par(mar=c(3,3,3,3))
  pie(mypie$Percent,labels = NA,col = mypie$cols, 
      main=paste("Cluster ",i,"\n(n=",sum(mypie$Freq),")", sep=""),
      border = 0,radius = 1)
  
  piepath=paste("output_files/revised_MAG/ward_cluster_k3_pie",i,".svg",sep="")
  mypie <- ggplot(mypie, aes(x="population",y=Freq, fill=Pop)) + 
    geom_col(colour="black")+ 
    geom_text(aes(label=as.character(Freq)), size=7, position=position_stack(vjust=0.5))+
    coord_polar(theta = "y") + 
    scale_fill_manual(values=mycols) + 
    theme_void(base_size = 20) + 
    ggtitle(paste("Cluster",i))
  mypie
  ggsave(piepath, plot=mypie,device = "svg", scale = 5,
         width = 5,height = 5,units = "cm",dpi = 300)

  }

# percentages
tab <- ftable(cluster.df$population,cluster.df$clust_k3)
round(prop.table(tab,margin=1)*100,0) -> cluster_countprop
tab <- ftable(cluster.df$urbanisation,cluster.df$clust_k3)
round(prop.table(tab,margin=1)*100,0) -> cluster_countLs

sink("output_files/revised_MAG/ward_cluster_k3_popcounts.txt")
print(tab)
print(cluster_countprop)
print(cluster_countLs)
sink()


