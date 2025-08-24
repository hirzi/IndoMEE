rm(list=ls())

setwd("C:/Users/caf77/OneDrive - University of Cambridge/Documents/Analysis/MOBILE/0_Pilot_reviewed/")

# load libraries
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(reshape2)
library(dplyr)
library(vegan)
source("pretty_plot_function_revised.R")

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
df <- read.delim("input_files/abundance_tables_raw/original/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv", header = T)
mypop <- readRDS("input_files/mypop_wreplicates.rds")
sampleid <- unique(as.character(mypop$sampleid))
pop.ord <- levels(mypop$population)

# Principal Coordinate Analysis
# With SPMP
# input comm
df <- read.delim("input_files/abundance_tables_raw/original/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv", header = T)
sampleid <- readRDS("input_files/sampleid_n116.rds")
sporeids <- readRDS("input_files/sampleid_spore.rds")
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
mycols <- read.delim("input_files/MAGs_pop_mycols.tsv", header=T)
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

# preview
pretty_PCA(cap, axes=c(1,2), metadata = env.df, 
           colour.order = pop.ord, colour.by = "population", colour.scheme = mycols,
           shape.by = "urbanisation", shape.order=names(myshape),shape.scheme = myshape)

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


#################
### Figure 2B ###
#################
# save plot
p <- pretty_PCA(cap, axes=c(1,2), metadata = env.df, 
                colour.order = pop.ord, colour.by = "population", colour.scheme = mycols,
                shape.by = "urbanisation", shape.order = names(myshape),shape.scheme = myshape)

saveRDS(p, "figout/Figure2/Fig2B.rds")
p <- readRDS("figout/Figure2/Fig2B.rds")
plot(p)

ggsave(plot=p,
       filename = "figout/Figure2/Fig2B.png", 
       device = "png", scale=4, dpi = 600, width = 5, height = 5, units = "cm")

ggsave(plot=p,
       filename = "figout/Figure2/Fig2B.svg", 
       device = "svg", scale=5, dpi = 600, width = 5, height = 5, units = "cm")

# save data
saveRDS(comm.df, "output_files/revised_MAG/pmv.comm.MAGs.indospmp.rds")
samples <- rownames(reads)
saveRDS(samples, "output_files/revised_MAG/pmv.sample.indospmp.rds")
saveRDS(env.df, "output_files/revised_MAG/pmv.metadata.indospmp.rds")

#########################
### IndoMEE only PCoA ###
#########################
# input comm
df <- read.delim("input_files/abundance_tables_raw/original/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv", header = T)
reads <- df[,colnames(df)%in%sampleid]
rownames(reads) <- df$Genome
reads <- convert_to_CPM(reads)
reads <- prevalence_filter(reads,2)
reads <- as.matrix(t(reads))

# input metadata
env.df <- readRDS("input_files/mypop_wSPMP.rds")
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
pop.ord.id <- pop.ord[pop.ord %in% unique(env.df$population)]
mycols <- read.delim("input_files/MAGs_pop_mycols.tsv", header=T)
mycols <- mycols[mycols$group %in% pop.ord.id,]
col.ord <- mycols$group
mycols <- mycols$mycols
names(mycols) <- col.ord

# preview
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

##################
### Figure S2B ###
##################
p <- pretty_PCA(cap, axes=c(1,2), metadata = env.df, 
                colour.order = pop.ord, colour.by = "population", colour.scheme = mycols,
                shape.by = "urbanisation", shape.order = names(myshape),shape.scheme = myshape)

saveRDS(p, "figout/FigureS2/FigS2-IndomeePCoA.rds")
p <- readRDS("figout/FigureS2/FigS2-IndomeePCoA.rds")
plot(p)


ggsave(plot=p,
       filename = "figout/FigureS2/FigS2-IndomeePCoA.png", 
       device = "png", scale=3, dpi = 1200, width = 5, height = 5, units = "cm")

ggsave(plot=p,
       filename = "figout/FigureS2/FigS2-IndomeePCoA.svg", 
       device = "svg", scale=3, dpi = 1200, width = 5, height = 5, units = "cm")



df <- read.delim("input_files/abundance_tables_raw/original/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv", header = T)
mypop <- readRDS("input_files/mypop_wreplicates.rds")
sampleid <- unique(as.character(mypop$sampleid))
pop.ord <- levels(mypop$population)


# save data
saveRDS(comm.df, "output_files/revised_MAG/pmv.comm.MAGs.indomee.rds")
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

# IndoMEE + SPMP
comm.df <- readRDS("output_files/revised_MAG/pmv.comm.MAGs.indospmp.rds")
comm.dist <- vegdist(comm.df+1, method="aitchison")
env.df <- readRDS("output_files/revised_MAG/pmv.metadata.indospmp.rds")
cap <- capscale(comm.dist~1, data=env.df, comm = comm.df)
centroids_IndoSPMP <- find_centroids(cap, env.df$population)

# IndoMEE
comm.df <- readRDS("output_files/revised_MAG/pmv.comm.MAGs.indomee.rds")
comm.dist <- vegdist(comm.df+1, method="aitchison")
env.df <- readRDS("output_files/revised_MAG/pmv.metadata.indomee.rds")
cap <- capscale(comm.dist~1, data=env.df, comm = comm.df)

centroids_Indomee <- find_centroids(cap, env.df$population)



##########################
### Sample Clustering ####
##########################
library(vegan)
library(ape)
library(clusterSim)
library(mclust)
library(ggplot2)
library(ggdendro)
library(ggpubr)

## Load Full data
comm.df <- readRDS("output_files/revised_MAG/pmv.comm.MAGs.indospmp.rds")
comm.dist <- vegdist(comm.df+1, method="aitchison")
mypop <- readRDS("output_files/revised_MAG/pmv.metadata.indospmp_revised.rds")
mypop$urbanisation <- factor(mypop$cLifestyle, c("remote","rural","urban"))

#1# Initial Clustering (IndoMEE + SPMP)
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
ward.clust.plot.SPMP <- ward.clust.plot
ggsave(plot = ward.clust.plot,"figout/additional_figs/phylo_SIN.png",device = "png",dpi = 600,width = 15, height = 2,units = "cm",bg = NULL, scale = 2.5)
saveRDS(ward.clust.plot, "figout/additional_figs/phylo_SIN.rds")

myclust <- as.phylo(mycluster) # generate tree file
write.tree(myclust,file = "figout/additional_figs/IndoMEE+SPMP.trefile",digits = 5)

# Optimise k
#2A# Elbow method to set tentative k
# Determine a tentatibe kmax value (k = number of cluster)
k.max=15
wss <- sapply(1:k.max, 
              function(k){kmeans(comm.dist, k, nstart=50,iter.max = 15 )$tot.withinss})
plot(1:k.max, wss, 
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares",
     main="Elbow Method, tentative point at k=4")
abline(v = 4, lty =2) # tentative elbow point, based on the rate of reduction of total-wss

#2B# Gap statistics # IndoMEE + SPMP samples
# Determine which k-value gives the most distinctive cluster (between k1 - kmax)
# To find the k with maximum distinction between clusters
# based on Tibshirani, Walther and Hastie gap index (https:/academic.oup.com/jrsssb/article/63/2/411/7083348)
# calculate gap between k=4, k=3, and k=2 (take the maximum gap)
clall <- cbind(cutree(mycluster,k=4), cutree(mycluster,k=3))
gap_k4_k3<- index.Gap(comm.df+1, clall,B=20,method = "ward.D2",centrotypes="centroids")
clall <- cbind(cutree(mycluster,k=3), cutree(mycluster,k=2))
gap_k3_k2 <- index.Gap(comm.df+1, clall,B=20,method = "ward.D2",centrotypes="centroids")
clall <- cbind(cutree(mycluster,k=2), 1)
gap_k2_k1 <- index.Gap(comm.df+1, clall,B=20,method = "ward.D2",centrotypes="centroids")
gap_stat <- c(gap_k4_k3$gap, gap_k3_k2$gap, gap_k2_k1$gap)
names(gap_stat) <- c("k4","k3","k2")

sink("figout/additional_figs/cluster_analysis_gap_stat.txt")
print("Gap Stat for all samples (IndoMEE + SPMP)")
gap_stat
gap_stat[which(gap_stat==max(gap_stat))]
sink()

#2C# Gap Statistics on only IndoMEE samples
# As any clustering methods are dataset-dependent, we want to make sure how well the clusters based on the IndoMEE+SPMP dataset performs on samples from IndoMEE-only
indomee.id <-  as.character(readRDS("output_files/revised_MAG/pmv.metadata.indomee_revised.rds")$sampleid)
comm.indomee <- comm.df[indomee.id,]
clall <- cbind(cutree(mycluster,k=4), cutree(mycluster,k=3))[indomee.id,]
gap_k4_k3<- index.Gap(comm.indomee+1, clall,B=20,method = "ward.D2",centrotypes="centroids")
clall <- cbind(cutree(mycluster,k=3), cutree(mycluster,k=2))[indomee.id,]
gap_k3_k2 <- index.Gap(comm.indomee+1, clall,B=20,method = "ward.D2",centrotypes="centroids")
clall <- cbind(cutree(mycluster,k=2), 1)[indomee.id,]
gap_k2_k1 <- index.Gap(comm.indomee+1, clall,B=20,method = "ward.D2",centrotypes="centroids")
gap_stat <- c(gap_k4_k3$gap, gap_k3_k2$gap, gap_k2_k1$gap)
names(gap_stat) <- c("k4","k3","k2")

sink("figout/additional_figs/cluster_analysis_gap_stat.txt",append = T)
print("Gap Stat for IndoMEE Only")
gap_stat
gap_stat[which(gap_stat==max(gap_stat))]
sink()

# save clusters as dataframe
clust_k4 <- cutree(mycluster,4)
clust_k3 <- cutree(mycluster,3)
clust_k2 <- cutree(mycluster,2)

cluster.df <- data.frame(mypop[mypop$sampleid%in%names(clust_k2),],
                         cbind(clust_k4,
                               clust_k3,
                               clust_k2))
cluster.df.SPMP<- cluster.df

# note: IndoMEE+SPMP dataset optimise at k=2 

#3# Compare with clusters generated on IndoMEE-only dataset
# load data
comm.df <- readRDS("output_files/revised_MAG/pmv.comm.MAGs.indomee.rds")
comm.dist <- vegdist(comm.df+1, method="aitchison")
mypop <- readRDS("output_files/revised_MAG/pmv.metadata.indomee_revised.rds")
mypop$urbanisation <- factor(mypop$cLifestyle, c("remote","rural","urban"))

#3A# Initial clustering (Ward.D2 method)
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
ward.clust.plot.IndoMEE <- ward.clust.plot
ggsave(plot = ward.clust.plot,"figout/additional_figs/phylo_IndoMEE.png",device = "png",dpi = 600,width = 15, height = 2,units = "cm",bg = NULL, scale = 2.5)
saveRDS(ward.clust.plot, "figout/additional_figs/phylo_IndoMEE.rds")

myclust <- as.phylo(mycluster) # generate tree file
write.tree(myclust,file = "figout/additional_figs/IndoMEE.trefile",digits = 5)

#3B# Elbow Method (determine Kmax)
k.max=15
wss <- sapply(1:k.max, 
              function(k){kmeans(comm.dist, k, nstart=50,iter.max = 15 )$tot.withinss})
plot(1:k.max, wss, 
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares",
     main="Elbow Method, tentative point at k=4")
abline(v = 4, lty =2) # tentative elbow point, based on the rate of reduction of total-wss

#3C# Gap Statistics 
clall <- cbind(cutree(mycluster,k=4), cutree(mycluster,k=3))
gap_k4_k3<- index.Gap(comm.df+1, clall,B=20,method = "ward.D2",centrotypes="centroids")
clall <- cbind(cutree(mycluster,k=3), cutree(mycluster,k=2))
gap_k3_k2 <- index.Gap(comm.df+1, clall,B=20,method = "ward.D2",centrotypes="centroids")
clall <- cbind(cutree(mycluster,k=2), 1)
gap_k2_k1 <- index.Gap(comm.df+1, clall,B=20,method = "ward.D2",centrotypes="centroids")
gap_stat <- c(gap_k4_k3$gap, gap_k3_k2$gap, gap_k2_k1$gap)
names(gap_stat) <- c("k4","k3","k2")
gap_stat
gap_stat[which(gap_stat==max(gap_stat))]

# note: IndoMEE-only dataset is optimise at k=2

# make dataframe
clust_k4 <- cutree(mycluster,4)
clust_k3 <- cutree(mycluster,3)
clust_k2 <- cutree(mycluster,2)

cluster.df <- data.frame(mypop[mypop$sampleid%in%names(clust_k2),],
                         cbind(clust_k4,
                               clust_k3,
                               clust_k2))
cluster.df.IndoMEE <- cluster.df

#4# Cluster Optimisation with ARI
# Aim: Investigate which k-point in the IndoMEE+SPMP dataset that is most similar to K=2 in IndoMEE-only dataset.
# Method: Adjusted Rand Index (ARI) measure the similarity of cluster configuration
# ARI value: 1 = perfect match, 0 = random, <0 = worse than random.

library(mclust)

indomee.id <- rownames(cluster.df.IndoMEE) # select only IndoMEE samples
clust_compare <- data.frame(cluster.df.SPMP[indomee.id,c("sampleid","clust_k2","clust_k3","clust_k4")],
                            clust_k2_indoMEE=cluster.df.IndoMEE$clust_k2)

clust_vs <- c("clust_k2","clust_k3","clust_k4")


adjustedRandIndex(clust_compare$clust_k2, clust_compare$clust_k2_indoMEE) # Calculate Adjusted Rand Index (ARI)

ari_out <- lapply(clust_vs, function(x){
  clust1 <- clust_compare[,"clust_k2_indoMEE"]
  clust2 <- clust_compare[,x]
  ari <- adjustedRandIndex(clust1, clust2)
  return(paste("indoMEE_k2 vs", x, "; ARI =", ari))
})

print(ari_out)
max_ari <- max(as.numeric(lapply(strsplit(x = unlist(ari_out),split = " = "), function(x) x[[2]])))
ari_out[grep(max_ari, ari_out)]

sink("figout/additional_figs/cluster_analysis_gap_stat.txt", append = T)
print("Adjusted Rand Index (ARI) to evaluate cluster structure")
print(ari_out)
print("cluster most congruent (max):")
print(ari_out[grep(max_ari, ari_out)])
sink()

# Note: based on ARI, we should choose k=3 for the SPMP+IndoMEE dataset

#5# Plot Sample distribution
cluster.df <- cluster.df.SPMP
p <- ggplot(cluster.df, aes(x=population, y=clust_k3, fill=population)) + theme_classic2(base_size = 12) +
  theme(panel.border = element_rect(fill = NA))+
  geom_vline(xintercept = seq(1.5,10.5,1))+
  geom_hline(yintercept = c(1.5,2.5,3.5))+
  geom_point(size=4, stat = "identity", pch=21, position = position_jitter(width = 0.45)) +
  scale_fill_manual(values = mycols)+
  scale_y_continuous(expand = c(0,0))
p
ggsave(p, filename = "figout/additional_figs/cluster_vs_urbanisation.svg",device = "png", height = 5, width = 5, units = "cm",dpi = 600)

#k3
tab <- ftable(cluster.df$population,cluster.df$clust_k3)
tab
chisq.test(tab,correct = T,rescale.p = T) #,B = 99999,simulate.p.value = T)

round(100*prop.table(tab,margin = 1),0) -> pop_clust
pop_clust # note that 2 = cluster3 ; 3 = cluster2

tab <- ftable(cluster.df$cLifestyle,cluster.df$clust_k3)
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
  
  piepath=paste("figout/additional_figs/ward_cluster_k3_pie",i,".svg",sep="")
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

sink("figout/additional_figs/ward_cluster_k3_popcounts.txt")
print(tab)
print(cluster_countprop)
print(cluster_countLs)
sink()

# Redraw tree with cluster
library(dplyr)
comm.df <- readRDS("output_files/revised_MAG/pmv.comm.MAGs.indospmp.rds")
comm.dist <- vegdist(comm.df+1, method="aitchison")
mypop <- readRDS("output_files/revised_MAG/pmv.metadata.indospmp_revised.rds")
mypop$urbanisation <- factor(mypop$cLifestyle, c("remote","rural","urban"))

mycluster <- hclust(comm.dist,method = "ward.D2")
clust <- cutree(mycluster, k=3)
cluster_df <- data.frame(label = names(clust), cluster = factor(clust))


hcdata <- dendro_data(mycluster, type="rectangle")
hcdata$labels <- hcdata$labels %>%
  left_join(cluster_df, by = "label")
mypop$label <- mypop$sampleid
hcdata$labels <- hcdata$labels %>%
  left_join(mypop[,c("label","population","urbanisation")], by="label")

tipsdata <- label(hcdata)
tipsdata$sampleid <- tipsdata$label

centers <- hcdata$labels %>%
  group_by(cluster) %>%
  summarize(x = min(x)+median(1:length(x)), y = 0)
centers$y <- centers$y + c(220,310,220)

ward.clust.plot <- ggplot() + 
  theme_void() +
  geom_segment(data=segment(hcdata), aes(x, y, xend=xend, yend=yend))+
  geom_point(data = tipsdata, 
             aes(x = x, y = y,col=population, fill=population), size=1.4, pch=22) +
  scale_colour_manual(values = mycols)+
  geom_label(data = centers,
             aes(x = x, y = y, label = paste("Cluster", cluster)),
             size = 3, fill = "white", label.size = 0.2,hjust = 0.5) +
  scale_fill_manual(values=mycols)
ward.clust.plot


# save cluster data
saveRDS(mycluster, "figout/additional_figs/ward_clusters_SIN.rds")
saveRDS(ward.clust.plot, "figout/Figure2/Fig2C_phylo_SIN.rds")

tmp <- cluster.df.IndoMEE
pf <- grep("clust",colnames(tmp))
colnames(tmp)[pf] <- paste0("IndoMEE_",colnames(tmp)[pf])
cluster.df <- merge.data.frame(cluster.df.SPMP, tmp[,c(1,pf)], by="sampleid")
write.table(cluster.df, "figout/additional_figs/ward_clusters.tsv", row.names = F, quote = F, sep = "\t")
