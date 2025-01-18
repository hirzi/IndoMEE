# Run functionalEnrichmentAnalysis.R till generation of feat.all
# Try also for presence-absence genes df

# Load libraries
library(tidyverse)
library(plotly)
library(factoextra)
library(compositions)

# Define run parameters
# Define whether to collapse count data to binary i.e. presence-absence
binary <- FALSE
clr_transform <- TRUE
scale_PCA <- FALSE
pops <- "indo_sg" # indo or indo_sg
if(binary == FALSE) {text_binary <- "counts"} else if (binary == TRUE) {text_binary <- "binary"}
if(clr_transform == TRUE) {text_CLR <- "CLR_"} else if(clr_transform == TRUE) {text_CLR <- ""}
if(scale_PCA == FALSE) {text_scale <- "notScaled"} else if (scale_PCA == TRUE) {text_scale <- "scaled"}
if(pops == "indo") {text_pops <- ""} else if (pops == "indo_sg") {text_pops <- "_wSPMP"} 

# Read in and format relative abundance table
if (pops == "indo") {
  relAbund_raw <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assemblies/On_hostDecontaminated_fastqs/HybridAssembly_Comp50Cont5_draft2/summary/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv", header = TRUE, sep = "\t")
} else if (pops == "indo_sg") {
  relAbund_raw <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/metamap/Assemblies/On_hostDecontaminated_fastqs/HybridAssembly_Comp50Cont5_draft2/summary_wSPMP/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv", header = TRUE, sep = "\t")
}

first_sample <- "ASMDAI005"
sample_idx <- match(first_sample, colnames(relAbund_raw))
relAbund <- relAbund_raw[,c(sample_idx:ncol(relAbund_raw))]
relAbund <- relAbund[,-c(grep("POSMOCKCTRL", colnames(relAbund)))]
rownames(relAbund) <- relAbund_raw$original_bin
relAbund$original_bin <- rownames(relAbund)

# Format functional annotations (gene counts) dataframe
feat.all2 <- feat.all
if(binary == TRUE) {
  feat.all2[feat.all2 > 1] <- 1
}
feat.all.renamed <- as.data.frame(t(feat.all2))
rownames(feat.all.renamed) <- gsub("treatment_", "", rownames(feat.all.renamed))
rownames(feat.all.renamed) <- gsub("control_", "", rownames(feat.all.renamed))
feat.all.renamed$original_bin <- rownames(feat.all.renamed)

# Merge relative abundance df with functional (gene counts) df
relAbund_feat <- merge(relAbund, feat.all.renamed, by = "original_bin")

# Remove unnecessary columns and samples
rownames(relAbund_feat) <- relAbund_feat$original_bin
relAbund_feat <- relAbund_feat[,-1]
first_sample <- "ASMDAI005"
last_sample <- "PBSHF003"
sampleFirst_idx <- match(first_sample, colnames(relAbund_feat))
sampleLast_idx <- match(last_sample, colnames(relAbund_feat))

# Get abundance-weighted gene counts across samples
for (i in sampleFirst_idx:sampleLast_idx) {
  relAbund_sample <- relAbund_feat[,c(i,(sampleLast_idx+1):ncol(relAbund_feat))]
  colnames(relAbund_sample)[1] <- "abundance"
  relAbund_sample_weighted <- relAbund_sample %>% 
    mutate(
      across(everything(), ~ .x* abundance)
    )
  relAbund_sample_weighted <- relAbund_sample_weighted[, -1]
  relAbund_sample_weightedSums <- as.data.frame(colSums(relAbund_sample_weighted))
  colnames(relAbund_sample_weightedSums) <- colnames(relAbund_feat)[i]
  if(i == 1) {
    relAbund_samples_weightedSums <- relAbund_sample_weightedSums
  } else if (i  > 1) {
    relAbund_samples_weightedSums <- cbind(relAbund_samples_weightedSums, relAbund_sample_weightedSums)
  }
}

relAbund_samples_weightedSums_t <- as.data.frame(t(relAbund_samples_weightedSums))
# CLR tranform abundances (think whether to do this or not)
if(clr_transform == TRUE) {
  #rowSums(relAbund_transposed) # should be ~1 i.e. matrix or data frame where each row is a composition
  # log(0) is undefined (= -Inf). Note that compositions::clr return 0 for log(0) through if/else statement, however, we want to avoid this behavior (since it is incorrect). 
  # To work around, we add an offset to the data (here, an order of magnitude smaller than the dataset's minimum value). For reference, see: http://mixomics.org/mixmc/mixmc-preprocessing/
  if(!is.na(match(0,relAbund_samples_weightedSums_t))) {
    min_val <- 10^(floor(log10(min(relAbund_samples_weightedSums_t[relAbund_samples_weightedSums_t > 0]))))
    relAbund_samples_weightedSums_t <- relAbund_samples_weightedSums_t + (min_val/10)
  }
  relAbund_samples_weightedSums_t <- as.data.frame(clr(relAbund_samples_weightedSums_t))
}

# Calculate PCA
if(scale_PCA == TRUE) {
  pc <- prcomp(relAbund_samples_weightedSums_t, center = TRUE, scale. = TRUE)
} else if (scale_PCA == FALSE) {
  pc <- prcomp(relAbund_samples_weightedSums_t, center = TRUE, scale. = FALSE)
}

# Check out eigenvalues
# pdf(file = paste0("/Users/hl636/Desktop/", feat, "_", text_binary, "_", text_CLR, text_scale, text_pops, "_eigenValues.pdf"), width = 10, height = 5)
# fviz_eig(pc)
# dev.off()
pc_summary <- summary(pc)
PC1_var <- pc_summary$importance[2,][1]*100
PC2_var <- pc_summary$importance[2,][2]*100

# Format PC dataframe
pc_df <- as.data.frame(pc$x)
pc_df$pop <- substring(rownames(pc_df), 1,5)
pc_df$pop[grep("BA0", pc_df$pop)] <- "Bali Denpasar"
pc_df$pop <- gsub("ASMDA", "Asmat", pc_df$pop)
pc_df$pop <- gsub("BAPDW", "Bali Pedawa", pc_df$pop)
pc_df$pop <- gsub("BRUSB", "Basap", pc_df$pop)
pc_df$pop <- gsub("MALHU", "Punan Tubu Hulu", pc_df$pop)
pc_df$pop <- gsub("MALRP", "Punan Tubu Respen", pc_df$pop)
pc_df$pop <- gsub("MALSU", "Punan Aput", pc_df$pop)
pc_df$pop <- gsub("PBSHF", "Punan Batu", pc_df$pop)
# Add Singapore pops
if (pops == "indo_sg") {
  chinese <- c("ERR7671879", "ERR7671880", "ERR7671888", "ERR7671915", "ERR7671920", "ERR7671922", "ERR7671925", "ERR7671933", "ERR7671935", "ERR7671938", "ERR7671941", "ERR7671950", "ERR7671954", "ERR7671964", "ERR7671975", "ERR7671877", "ERR7671883", "ERR7671895", "ERR7671899", "ERR7671901", "ERR7671902", "ERR7671912", "ERR7671913", "ERR7671918", "ERR7671921", "ERR7671923", "ERR7671942", "ERR7671960", "ERR7671961", "ERR7671965", "ERR7671967", "ERR7671981", "ERR7671878", "ERR7671881", "ERR7671884", "ERR7671892", "ERR7671894", "ERR7671898", "ERR7671910", "ERR7671919", "ERR7671932", "ERR7671934", "ERR7671940", "ERR7671945", "ERR7671948", "ERR7671949", "ERR7671952", "ERR7671953", "ERR7671956", "ERR7671957", "ERR7671959", "ERR7671973", "ERR7671980")
  indian <- c("ERR7671874", "ERR7671875", "ERR7671886", "ERR7671889", "ERR7671907", "ERR7671930", "ERR7671946", "ERR7671976", "ERR7671979", "ERR7671982", "ERR7671887", "ERR7671891", "ERR7671896", "ERR7671903", "ERR7671904", "ERR7671905", "ERR7671906", "ERR7671909", "ERR7671924", "ERR7671927", "ERR7671951", "ERR7671958", "ERR7671962", "ERR7671911", "ERR7671943", "ERR7671955", "ERR7671963", "ERR7671966", "ERR7671968", "ERR7671974")
  malay <- c("ERR7671882", "ERR7671885", "ERR7671890", "ERR7671893", "ERR7671908", "ERR7671917", "ERR7671928", "ERR7671939", "ERR7671947", "ERR7671971", "ERR7671978", "ERR7671876", "ERR7671897", "ERR7671914", "ERR7671926", "ERR7671929", "ERR7671931", "ERR7671944", "ERR7671969", "ERR7671972", "ERR7671977", "ERR7671900", "ERR7671916", "ERR7671936", "ERR7671937", "ERR7671970")
  Sg_pop_names <- c("Singapore - Chinese", "Singapore - Indian", "Singapore - Malay")
  Sg_idx <- list(chinese, indian, malay)
  for(j in seq(1,3)) {
    idx_list <- match(sort(Sg_idx[[j]]), rownames(pc_df))
    for(k in idx_list) {
      pc_df$pop[k] <- Sg_pop_names[j]
    }
  }
}
pc_df$pop <- as.factor(pc_df$pop)

# Plot PCA
pop_cols <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Pilot ms/pop_colours.tsv", header = TRUE, sep = "\t", comment.char="") 
if (pops == "indo") {
  col_pops <- c("Asmat" = "#FF6347", "Bali Denpasar" = "#2E8B57", "Bali Pedawa" = "#54FF9F", "Basap" = "#FFB90F", "Punan Tubu Hulu" = "#904E9F", "Punan Tubu Respen" = "#DF96B4", "Punan Aput" = "#B494D5", "Punan Batu" = "#5B3794")
} else if (pops == "indo_sg") {
  col_pops <- c("Asmat" = "#FF6347", "Bali Denpasar" = "#2E8B57", "Bali Pedawa" = "#54FF9F", "Basap" = "#FFB90F", "Punan Tubu Hulu" = "#904E9F", "Punan Tubu Respen" = "#DF96B4", "Punan Aput" = "#B494D5", "Punan Batu" = "#5B3794", "Singapore - Indian" = "#27408B", "Singapore - Malay" = "#4876FF", "Singapore - Chinese" = "#A6CCFF")
}
#pdf(file = paste0("/Users/hl636/Desktop/", feat, "_", text_binary, "_", text_CLR, text_scale, text_pops, ".pdf"), width = 10, height = 10*PC2_var/PC1_var*0.75) # scaling aspect ratio proportional to PC1:PC2 variance explained, with additional coeff to account for legend
# pdf(file = paste0("/Users/hl636/Desktop/", feat, "_", text_binary, "_", text_CLR, text_scale, text_pops, ".pdf"), width = 10, height = 5)
pl <- ggplot() + geom_point(data=pc_df, size = 3, aes_string(x="PC1", y="PC2", color=pc_df$pop)) + theme_bw() +  scale_color_manual(name = "POPULATION", values = col_pops) + labs(x = paste0("PC1 (", PC1_var,"%)"), y = paste0("PC2 (", PC2_var,"%)"))
#print(pl)
# dev.off()

# Make an interactive plotly plot
#pl_plotly <- ggplotly(pl)
#pl_plotly


##### For nicer plot #####
library(ggplot2)
library(ggpubr)
library(gridExtra)

# Read in data
setwd("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Genome-based reference and de-novo assemblies/Functional enrichment/")
mypop <- readRDS("mypop.rds")

# limit settings normal
#xlims <- c(-6,4); ylims <- c(-5,5)
xlims <- c(min(pc_df$PC1)*1.2,max(pc_df$PC1)*1.2); ylims <- c(min(pc_df$PC2)*1.2,max(pc_df$PC2)*1.2)

# xylimits settings presence/absence
# xlims <- c(-2,2); ylims <- c(-2,2)

# select PC
PC <- c(1,2)

# prepare plot data for collation
xycoord <- cbind(pc_df[, c(1,2)], rownames(pc_df))
rownames(xycoord) <- NULL
colnames(xycoord) <- c("Xaxis", "Yaxis", "sampleid")
pco.coord <- xycoord # save for later
var <- pc_summary$importance[2,]

# format axis labels
axis_labs <- paste("PC",1:length(names(var)), " (", round(var*100,1),"%)", sep="")

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


# #### Look at PC loadings ####
# pc_loadings <- as.data.frame(pc$rotation)
# par(mfrow=c(1,1))
# hist(pc_loadings$PC1, 100)
# 
# # Grep sporulation indices in pc loadings df
# sp_1_idx <- grep("spore", rownames(pc_loadings))
# sp_2_idx <- grep("sporu", rownames(pc_loadings))
# sp_idx <- sort(unique(c(sp_1_idx,sp_2_idx)))
# 
# # Add lines for (all) sporulations genes
# sp_loadings_all <- pc_loadings$PC1[sp_idx]
# abline(v=sp_loadings_all, col="blue")
# 
# # Add lines for sporulation genes enriched in novel MAGs
# # Read in MAG metadata
# MAG.metadata2 <- read.table("/Users/hl636/Documents/Hirzi/Cambridge/Metagenomics/Phylogenetics/HybridAssembly_Comp50Cont5_draft2/Metadata_HybridAssembly_Comp50Cont5_draft2_revised25Nov_MaaslinNeedsToBeUpdated.txt", comment.char="", sep = "\t", header = TRUE)
# # Get MAG presence-absence of significantly associated sporulation-associated features
# feat.all.temp <- feat.all
# feat.all.temp$feature <- rownames(feat.all.temp)
# feat.sel <- merge(res.sig, feat.all.temp, by="feature")
# feat.sel$feature <- str_split_fixed(feat.sel$feature, ";", n = Inf)[,1]
# feat.sel1 <- feat.sel[grep("spore", feat.sel$description), ]
# feat.sel2 <- feat.sel[grep("sporu", feat.sel$description), ]
# feat.sel <- merge(feat.sel1, feat.sel2, all=TRUE)
# # Add lines
# sp_sel_idx <- match(feat.sel$description, rownames(pc_loadings))
# sp_loadings_sel <- pc_loadings$PC1[sp_sel_idx]
# abline(v=sp_loadings_sel, col="red")
# 
# # Plot histogram of MAG and KEGG ortholog relative abundances
# relAbund2 <- relAbund_raw[,c(sample_idx:ncol(relAbund_raw))]
# relAbund2 <- relAbund2[,-c(grep("POSMOCKCTRL", colnames(relAbund2)))]
# rownames(relAbund2) <- relAbund_raw$original_bin
# total_relative_abundance_MAG <- rowSums(relAbund2)/(sum(rowSums(relAbund2)))
# total_relative_abundance_KEGG <- rowSums(relAbund_samples_weightedSums)/(sum(rowSums(relAbund_samples_weightedSums)))
# par(mfrow=c(2,1))
# hist(total_relative_abundance_MAG, breaks = 400)
# hist(total_relative_abundance_KEGG, breaks = 400)
