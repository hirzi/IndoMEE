
# Load libraries
library(vegan)
library(compositions)

# Set working directory
setwd("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/Functional enrichment/")

# Define parameters
pops <- "indo_sg" # indo or indo_sg
distance_metric <- "jaccard" # "jaccard", "aitchison", "robust.aitchison" or any vegdist method
ord_method <- "aitchison_capscale" # "aitchison_capscale" or "clr_pca"
if(ord_method == "clr_pca") {
  scale_PCA <- FALSE
}

# Read in data
mypop <- readRDS("mypop.rds")
sampleid <- readRDS("sampleid_final.rds")
if(pops == "indo_sg") {
  df <- read.delim("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assemblies/On_hostDecontaminated_fastqs/HybridAssembly_Comp50Cont5_draft2/summary_wSPMP/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv")
} else if(pops == "indo") {
  df <- read.delim("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assemblies/On_hostDecontaminated_fastqs/HybridAssembly_Comp50Cont5_draft2/summary/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv")
  mypop <- mypop[mypop$population != "Malay_SIN",]
  mypop <- mypop[mypop$population != "Chinese_SIN",]
  mypop <- mypop[mypop$population != "Indian_SIN",]
  SPMP_samples_idx <- grep("ERR", sampleid)
  sampleid <- sampleid[-SPMP_samples_idx]
}

# Define run parameters
filter_MAGs <- TRUE
if(filter_MAGs==TRUE) {
  MAG_type <- "non_sporers" # sporers or non_sporers
  MAG_type_thres <- "counts" # counts or binary
  label <- paste0(MAG_type, "_", MAG_type_thres)
} else if (filter_MAGs==FALSE) {
  label <- "all_MAGs"
}

# Filter by MAGs (e.g. sporing / non-sporing)
df_filt <- df
if(filter_MAGs==TRUE) {
  if(MAG_type_thres=="binary"){
    # Get spo0A MAGs. First, run functionalEnrichmentAnalysis.R
    feat.spo0A <- feat.all[grep("spo0A", rownames(feat.all)), ]
    colnames(feat.spo0A) <- gsub("treatment_", "", colnames(feat.spo0A))
    colnames(feat.spo0A) <- gsub("control_", "", colnames(feat.spo0A))
    feat.spo0A <- as.data.frame(t(feat.spo0A))
    colnames(feat.spo0A) <- "spo0A"
    feat.spo0A$original_bin <- rownames(feat.spo0A)
    spore_MAGs <- feat.spo0A[feat.spo0A$spo0A==1,]$original_bin
    nonSpore_MAGs <- feat.spo0A[feat.spo0A$spo0A==0,]$original_bin
    spore_MAGs <- sample(spore_MAGs, min(length(spore_MAGs), length(nonSpore_MAGs)))
    nonSpore_MAGs <- sample(nonSpore_MAGs, min(length(spore_MAGs), length(nonSpore_MAGs)))
  } else if(MAG_type_thres=="counts"){
    # Alternatively, partition by spore counts. Get spore counts for Bacilotta A. First, run functionalEnrichmentAnalysis.R, followed by Sporulation_features_presenceAbsence_logisticRegression.R (till the appropriate row).
    # Consider parition e.g. top 200 vs bottom 200.
    p <- phyla[3]
    print(p)
    phylum_df <- phyla_df[phyla_df$phylum == p, ]
    feat.sel.prev.phylum <- merge(phylum_df, feat.sel.prev, by = "original_bin")
    feat.sel.sum <- feat.sel.prev.phylum[,c(2:(ncol(feat.sel.prev.phylum)-1))]
    feat.sel.sum <- feat.sel.sum[,-c(1)]
    feat.sel.sum$sum <- rowSums(feat.sel.sum)
    feat.sel.sum$prevalence <-feat.sel.prev.phylum$prevalence
    feat.sum.BacA <- as.data.frame(feat.sel.sum$sum)
    rownames(feat.sum.BacA) <- feat.sel.prev.phylum$original_bin
    colnames(feat.sum.BacA) <- c("sum")
    feat.sum.BacA$original_bin <- rownames(feat.sum.BacA)
    spore_MAGs <- feat.sum.BacA[feat.sum.BacA$sum>23,]$original_bin
    nonSpore_MAGs <- feat.sum.BacA[feat.sum.BacA$sum<10,]$original_bin
    spore_MAGs <- sample(spore_MAGs, min(length(spore_MAGs), length(nonSpore_MAGs)))
    nonSpore_MAGs <- sample(nonSpore_MAGs, min(length(spore_MAGs), length(nonSpore_MAGs)))
  }
  # And filter
  rownames(df_filt) <- df_filt$original_bin
  if (MAG_type == "sporers") {
    df_filt <- df_filt[rownames(df_filt) %in% spore_MAGs, ]
  } else if (MAG_type == "non_sporers") {
    df_filt <- df_filt[rownames(df_filt) %in% nonSpore_MAGs, ]
  }
}


# Filter by sample
comm.df <- df_filt[,colnames(df_filt) %in% sampleid]
comm.df <- comm.df*10^6
rownames(comm.df) <- df_filt$Genome
comm.df <- data.frame(t(comm.df))
comm.df <- comm.df[sampleid,]
comm.df <- comm.df[rowSums(comm.df)>0,]
table(colSums(comm.df)==0)
dim(comm.df)
pc.sampleid <- rownames(comm.df)

# metadata
env.df <- mypop[mypop$sampleid %in% sampleid,]
removed.samples <- mypop$sampleid[!mypop$sampleid %in% sampleid]
rownames(env.df) <- env.df$sampleid
env.df <- env.df[,c("sampleid","population")]
env.df <- env.df[sampleid,]
data.frame(table(env.df$population))

############
### PCoA ###
############

if(ord_method == "clr_pca") {
  rowSums(comm.df) # should be ~1x10ˆ6
  # log(0) is undefined (= -Inf). Note that compositions::clr return 0 for log(0) through if/else statement, however, we want to avoid this behavior (since it is incorrect). 
  # To work around, we add an offset to the data (here, an order of magnitude smaller than the dataset's minimum value), or 1 for counts data. For reference, see: http://mixomics.org/mixmc/mixmc-preprocessing/
  comm.df.clr <- as.data.frame(t(clr(comm.df + 1)))
  comm.df.clr.t <- as.data.frame(t(comm.df.clr))
  if(scale_PCA == TRUE) {
    pc <- prcomp(comm.df.clr.t, center = TRUE, scale. = TRUE)
  } else if (scale_PCA == FALSE) {
    pc <- prcomp(comm.df.clr.t, center = TRUE, scale. = FALSE)
  }
  pc_summary <- summary(pc)
  pc_df <- as.data.frame(pc$x)
  var <- pc_summary$importance[2,]
} else if(ord_method == "aitchison_capscale") {
  comm.dist <- vegdist(comm.df+1, method=distance_metric)
  cap <- capscale(comm.dist~1, comm=comm.df)
  eig <- cap$CA$eig
  var <- eig/sum(eig)
}

# presence absence filter
# raup.df <- comm.df
# raup.df <- raup.df[rowSums(raup.df) > 0,]
# raup.df[raup.df>0] <- 1
# comm.dist <- vegdist(raup.df, method="raup")
# pc.sampleid <- rownames(raup.df)
# dim(raup.df)

####################
# Draw pretty plot #
####################

library(ggplot2)
library(ggpubr)
library(gridExtra)

# set colours
mycols <- read.delim("MAGs_pop_mycols.tsv", header=T)

# select PC
PC <- c(1,2)

# prepare plot data for collation
if(ord_method == "clr_pca") {
  xycoord <- cbind(pc_df[, c(1,2)], rownames(pc_df))
  axis_labs <- paste("PC",1:length(names(var)), " (", round(var*100,1),"%)", sep="")
} else if(ord_method == "aitchison_capscale") {
  xycoord <- summary(cap)
  xycoord <- data.frame(xycoord$sites[,PC])
  xycoord$sampleid <- rownames(xycoord)
  axis_labs <- paste("PCo",1:length(names(var)), " (", round(var*100,1),"%)", sep="")
}
rownames(xycoord) <- NULL
pco.coord <- xycoord # save for later
colnames(xycoord) <- c("Xaxis","Yaxis","sampleid")

# limit settings normal
#lims <- c(-6,4); ylims <- c(-5,5)
xlims <- c(min(xycoord$Xaxis)*1.2,max(xycoord$Xaxis)*1.2); ylims <- c(min(xycoord$Yaxis)*1.2,max(xycoord$Yaxis)*1.2)
# xylimits settings presence/absence
# xlims <- c(-2,2); ylims <- c(-2,2)

# set population order
pop.ord <- c("Asmat_DAI","Punan_HUL","Punan_PBS","Punan_SUL","Punan_RPN",
             "Basap_BRU","Balinese_PDW","Balinese_DPS","Jakarta_JKT",
             "MOCK_CTR","NEG_CTR", "Dogs", "Oral",
             "Indian_SIN","Malay_SIN","Chinese_SIN")

# set environmental variables
ggdata <- merge.data.frame(xycoord, mypop, by="sampleid")

ggdata$population <- factor(ggdata$population, pop.ord)
ggdata <- ggdata[order(ggdata$population),]

ggdata$sampleid <- factor(ggdata$sampleid, unique(ggdata$sampleid))
rownames(ggdata) <-  NULL

# set axis title
xtitle <- axis_labs[PC[1]]
ytitle <- axis_labs[PC[2]]

# set color by country
mycols <- read.delim("MAGs_pop_mycols.tsv", header=T)
col.ord <- mycols$group
mycols <- mycols$mycols
names(mycols) <- col.ord
mycols <- c(mycols)
sincols <- c("#27408B","#4876FF","#CAE1FF")
names(sincols) <- c("Indian_SIN","Malay_SIN","Chinese_SIN")
mycols <- c(mycols,sincols)
mycols <- mycols[names(mycols)%in%pop.ord]

# make base plot
ggdata$population <- factor(ggdata$population, pop.ord)
a <- ggplot() +
  theme_bw() +
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept=0, lty="dashed")

# draw the rest of the plot
a <- a +
  geom_point(data=ggdata, aes(x=Xaxis, y=Yaxis, fill=population), col="black", size=3, pch=21) +
  geom_point(data=ggdata, aes(x=Xaxis, y=Yaxis, fill=population), col="black", size=3, pch=21) +
  scale_fill_manual(values=mycols) +
  scale_x_continuous(expand = c(0,0),
                     limits = xlims,
                     breaks = round(seq(xlims[1],xlims[2],length.out = 5)),
                     position = "top")+
  scale_y_continuous(expand = c(0,0),
                     limits = ylims,
                     breaks = round(seq(ylims[1],ylims[2],length.out = 5)),
                     position = "left") +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=12),
        plot.margin=margin(0,0,0,0,"cm"),
        legend.position = "none")+
  xlab(xtitle) + ylab(ytitle)


# PCo Group Boxplot X axis
pop_cols <- mycols[pop.ord] # set color by population

ggdata$population <- factor(ggdata$population, rev(pop.ord))
ggdata <- ggdata[order(ggdata$population),]
b <- ggplot(ggdata,
            aes(x=population, y=Xaxis, fill=population))+
  theme_bw()+
  geom_hline(yintercept = 0, lty="dashed") +
  geom_boxplot(width=0.8, outlier.size = 2, outlier.shape=21,
               outlier.fill = "white") +
  scale_fill_manual(values=pop_cols[pop.ord]) +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_blank(),
        axis.title.y = element_text(colour = "white",size=12),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin=margin(0,0,0,0,"cm"),
        legend.position="none")+
  scale_y_continuous(expand = c(0,0),
                     limits = xlims,
                     breaks = round(seq(xlims[1],xlims[2],length.out = 5)),
                     position = "left") +
  xlab("sample distance") + ylab(" ")+
  coord_flip()

# PCo Group Boxplot Y axis
ggdata$population <- factor(ggdata$population, pop.ord)
c <- ggplot(ggdata,
            aes(x=population, y=Yaxis, fill=population)) +
  theme_bw()+
  geom_hline(yintercept = 0, lty="dashed")+
  geom_boxplot(width=0.8, outlier.size = 2,
               outlier.shape=21, outlier.fill = "white") +
  scale_fill_manual(values=pop_cols) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x=element_text(colour = "white",size=12),
        plot.margin=margin(0,0,0,0,"cm"),
        legend.position = "none")+
  scale_y_continuous(expand = c(0,0),
                     limits = ylims,
                     breaks = round(seq(ylims[1],ylims[2],length.out = 5)),
                     position = "right")+
  scale_x_discrete(position="top")+
  xlab("sample distance") + ylab(" ")


# fix labels and extract legend

pop_count <- data.frame(table(mypop$population))
colnames(pop_count)  <- c("population","sample")
pop_count$population <- factor(pop_count$population, pop.ord)

poplab <-paste0(pop_count$population," (n=",pop_count$sample,")")
names(poplab) <- pop_count$population

ggdata$label <- c()
for(pop in names(poplab)){
  ggdata[ggdata$population == pop, "label"] <- poplab[pop]
}
ggdata <- ggdata[order(ggdata$population),]
ggdata$label <- factor(ggdata$label, unique(ggdata$label))

legend.col <-  mycols[as.character(unique(ggdata$population))]
names(legend.col) <- poplab[names(legend.col)]

p <- ggplot(ggdata, aes(x=label, y=Xaxis, fill=label))+ theme_bw() +
  geom_boxplot(width=0.8, outlier.size = 2, outlier.shape=21, outlier.fill = "white") +
  scale_fill_manual(values=legend.col) +
  theme(legend.text = element_text(size=10),
        legend.title = element_blank(),
        legend.key.size = unit(0.3,"cm"),
        legend.background = element_blank())
legend <- cowplot::get_legend(p)
legend <- as_ggplot(legend) + theme(plot.margin = margin(-1.5,1,0,0,"cm"))


# set plot layout
lay <- rbind(c(1,1,1,1,2,2),
             c(1,1,1,1,2,2),
             c(1,1,1,1,2,2),
             c(1,1,1,1,2,2),
             c(3,3,3,3,4,4),
             c(3,3,3,3,4,4))
# draw plot
print("Plot: Aitchison Distance for CLR-transformed Data")
grid.arrange(a,c,b,legend, layout_matrix=lay)


#### Look at PC loadings ####

# Read in novel MAG metadata
metadata_raw <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Phylogenetics/HybridAssembly_Comp50Cont5_draft2/Metadata_HybridAssembly_Comp50Cont5_draft2_revised25Nov_MaaslinNeedsToBeUpdated.txt", header=TRUE, comment.char = "", sep = "\t")
novelMAGs <- metadata_raw[metadata_raw$novelty==1,]$genome

# Plot histogram
pc_loadings <- as.data.frame(summary(cap)$species)
#write.table(pc_loadings, "/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/Functional enrichment/allPhyla_allMAGs_IndowSPMP_MAGAbundPCoA_Jaccard_loadings_revised.txt", row.names = TRUE, quote = FALSE, sep = "\t")
hist(pc_loadings$MDS1, 100)

# Grep novel MAGs' indices in pc loadings df
novelMAGs_idx <- match(novelMAGs, rownames(pc_loadings))
novelMAGs_loadings <- pc_loadings$MDS1[novelMAGs_idx]

# Add lines for (all) sporulations genes
abline(v=novelMAGs_loadings, col="blue")
