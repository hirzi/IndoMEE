library(ggplot2)
library(ggpubr)
library(ggtext)
library(gridExtra)

## custom function
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

select_taxon_by_coef <- function(maaslin.output, coef.column, n){
  maaslin.output <- maaslin.output[order(maaslin.output[, coef.column], decreasing = T),]
  select <- maaslin.output$feature[1:n]
  select}

###################
# Append taxonomy #
###################
full.taxonomy <- read.delim("output_files/revised_MAG/MAASLIN2/indomee_full_taxonomy.txt")
colnames(full.taxonomy) <- tolower(colnames(full.taxonomy))
full.taxonomy$species_old <- full.taxonomy$species
full.taxonomy$species <- full.taxonomy$species_deduplicated
full.taxonomy$feature <- full.taxonomy$original_bin

# full dataset
data <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/cont_p20/all_results.tsv", header=T)
data <- merge.data.frame(data, full.taxonomy[,c("feature","genome","original_bin","classification","phylum","family","genus","species")],by = "feature")
write.table(data, "output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/cont_p20/all_results_wtax.tsv",quote = F,row.names = F,sep = "\t")
data$feature <- data$genome
data$feature <- gsub("MGYG","INDOMEE", data$feature)
write.table(data, "output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/cont_p20/all_results_MGYG_wTax.tsv",quote = F,row.names = F,sep = "\t")

tmp <- data[data$metadata=="lifestyle",]
table(duplicated(tmp$feature))


# Indonesia-only dataset
data <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/cont_p20_indomee2/all_results.tsv", header=T)
data <- merge.data.frame(data, full.taxonomy[,c("feature","genome","original_bin","classification","phylum","family","genus","species")],by = "feature")
write.table(data, "output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/cont_p20_indomee2/all_results_wtax.tsv",quote = F,row.names = F,sep = "\t")
data$feature <- data$genome
data$feature <- gsub("MGYG","INDOMEE", data$feature)
write.table(data, "output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/cont_p20_indomee2/all_results_MGYG_wTax.tsv",quote = F,row.names = F,sep = "\t")

tmp <- data[data$metadata=="lifestyle",]
table(duplicated(tmp$feature))

# edit feature in the taxonomy key
full.taxonomy$feature <- gsub("MGYG","INDOMEE",full.taxonomy$genome)

#####################
# Analyse Full Data #
#####################
# read input
input.df <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/input.comm_p20.tsv", header=T)
metadata <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/input.metadata_p20.tsv", header=T)
dim(input.df)
df <- merge.data.frame(metadata, input.df, by="sampleid")
rownames(df) <- df$sampleid
RA <- input.df

data <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/cont_p20/all_results_MGYG_wTax.tsv")
data$slope <- sign(data$coef)

signif.df <- data[data$qval<0.05,]
signif.df <- signif.df[signif.df$value=="lifestyle",]
signif.df <- signif.df[order(signif.df$coef),]

length(unique(data$feature))
length(unique(signif.df$feature))
100*length(unique(signif.df$feature))/561
100*length(unique(signif.df$feature))/length(unique(data$feature))
# chisqr by signs
data[data$qval >0.05,"slope"] <- 0

table(data$phylum, data$slope)
chisq.test(table(data$phylum, data$slope))
chisq.test(table(data$family, data$slope))
chisq.test(table(data$genus, data$slope))


############
## Top 60 ##
############
heat.data <- data[order(data$coef, decreasing = T),]
heat.data <- heat.data[heat.data$metadata=="lifestyle",]
heat.data <- heat.data[heat.data$qval<0.05,]
heat.data$abs.coef <- abs(heat.data$coef)

# select top 60
select <- select_taxon_by_coef(heat.data,coef.column = "abs.coef",n = 60)
saveRDS(select, file="output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/cont_p20/top60_features.rds")

# Match with Indonesia-only dataset
df.indomee <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/cont_p20_indomee2/all_results_MGYG_wtax.tsv", header=T)
indomee.mags <- unique(df.indomee$feature)
indomee_extra_mags <- indomee.mags[!indomee.mags %in% data$feature]
length(indomee_extra_mags) # list of mags that only passed the p>0 filter in the Indonesian-only dataset

signif.df <- df.indomee[df.indomee$qval<0.05,]
signif.df <- signif.df[signif.df$value=="lifestyle",]
indomee_signif_mags <- unique(signif.df$feature)
length(unique(indomee_signif_mags))

spmp_signif_mags <- data[data$qval<0.05,]
spmp_signif_mags <- unique(spmp_signif_mags[spmp_signif_mags$value=="lifestyle","feature"])


spmp_signif_indoMatch_mags <-  spmp_signif_mags[spmp_signif_mags%in%indomee_signif_mags]
length(spmp_signif_indoMatch_mags)

100*length(spmp_signif_indoMatch_mags)/length(spmp_signif_mags)

library(ggvenn)
ggvenn(data=list(All_Samples=spmp_signif_mags,
                 IndoMEE=indomee_signif_mags))

# coeficient correlation and direction
coef.spmp <- data[data$metadata=="lifestyle" & (data$feature %in% spmp_signif_indoMatch_mags),c("feature","coef")]
coef.indomee <- df.indomee[df.indomee$metadata=="lifestyle" & (df.indomee$feature %in% spmp_signif_indoMatch_mags),c("feature","coef")]
colnames(coef.spmp)[2] <- "coef.spmp"
colnames(coef.indomee)[2] <- "coef.indomee"
coef.compare <- merge.data.frame(coef.spmp, coef.indomee)
coef.compare$same_direction <- ifelse(sign(coef.compare$coef.spmp) == sign(coef.compare$coef.indomee), T, F)
cor.test(coef.compare$coef.spmp, coef.compare$coef.indomee)
table(coef.compare$same_direction)

##########
# Bifido #
##########
x <- full.taxonomy$species[full.taxonomy$feature%in%select]
x[grep("Bifidobacterium",x)] # list in top 60

# Pearson corr
n <- c(full.taxonomy$original_bin[grep("Bifidobacterium",full.taxonomy$genus)],
       full.taxonomy$original_bin[grep("Treponema",full.taxonomy$genus)])
corr.mat <- data.frame(input.df[,colnames(input.df) %in% n])
colnames(corr.mat) <- unlist(lapply(colnames(corr.mat), function(x) {
  x <- full.taxonomy$species[which(full.taxonomy$original_bin==x)]
  return(x)}))
corr.mat <- corr.mat[,sort(colnames(corr.mat))]
corr.mat <- corr.mat[,-grep("dentium",colnames(corr.mat))]
psych::corr.test(corr.mat,method = "kendall") -> corr.out
round(corr.out$r[c(1:4),-c(1:4)],2) # rho
round(corr.out$p[c(1:4),-c(1:4)],4) # pval

# signif in Indonesia
x <- df.indomee[grep("Bifidobacterium",df.indomee$species),]
x <- x[x$metadata=="lifestyle",]
x$species[x$qval<0.05]

# abundance and occurance
indomee_sid <- readRDS("output_files/sampleid_n116.rds")
spmp_sid <- metadata$sampleid

Bifid <- full.taxonomy$original_bin[grep("Bifidobacterium", full.taxonomy$species)]
Bifid <- input.df[,c("sampleid",Bifid)]
Bifid$population <- NA
Bifid$population[grep("ERR",Bifid$sampleid)] <- "Singapore"
Bifid$population[-grep("ERR",Bifid$sampleid)] <- "Indonesia"
table(Bifid$population)

Bifid_compare <- apply(Bifid[,2:6],2, function(x) {wilcox.test(x ~ Bifid$population)})
names(Bifid_compare) <- unlist(lapply(names(Bifid_compare), function(x) {
  x <- full.taxonomy$species[which(full.taxonomy$original_bin==x)]
  return(x)}))
Bifid_compare

Bifid_summary <- do.call(rbind,lapply(2:6, function(x) {
  mag <- colnames(Bifid)[x]
  x <- Bifid[,x]
  prev <- x
  prev[prev>0] <- 1
  sin <- grep("ERR", Bifid$sampleid)
  prev <- c(Indonesia=sum(prev[-sin]), Singapore=sum(prev[sin]))
  i <- summary(x[-sin])
  s <- summary(x[sin])
  out <- rbind(i,s)
  si <- sum(x[-sin])
  ss <- sum(x[sin])
  out <- cbind(out, rbind(si,ss))
  colnames(out)[7] <- "sumTotal"
  
  spec <- full.taxonomy$species[full.taxonomy$original_bin==mag]
  out <- data.frame(prevalence=prev, out, species=spec, MAG=mag, population= c("Indonesia","Singapore"), row.names = NULL)
  return(out)
}))

Bifid_summary

for(i in unique(Bifid_summary$MAG)) {
  sp <- Bifid_summary[Bifid_summary$MAG==i&Bifid_summary$population=="Singapore","sumTotal"]
  ip <- Bifid_summary[Bifid_summary$MAG==i&Bifid_summary$population=="Indonesia","sumTotal"]
  rat <- ip / sp
  if(rat < 1) {rat <- -1 * (1/rat)}
  spc <- full.taxonomy$species[full.taxonomy$original_bin == i]
  out <- paste(spc, i, round(rat,3) , sep = " = ")
  out <- paste("Indonesia : Singapore fold change", out)
  print(out)
}

library(dplyr)
Bifid$population <- NULL
Bifid <- merge.data.frame(Bifid, metadata[,c("sampleid","population")], by="sampleid")
Bifid_summary2 <- do.call(rbind,lapply(2:6, function(x) {
  mag <- colnames(Bifid)[x]
  x <- Bifid[,x]
  group <- Bifid$population
  prev <- x
  prev[prev>0] <- 1
  data.frame(prev=prev,group=group) %>% group_by(group) %>% summarise(prevalence=sum(prev)) -> prev
  colnames(prev) <- c("population","prevalence")
  
  data.frame(relab=x,group=group) %>% group_by(group) %>% summarise(average=mean(relab,na.rm = T),
                                                                    median=median(relab, na.rm = T),
                                                                    minimum=min(relab, na.rm = T),
                                                                    maximum=max(relab, na.rm = T),
                                                                    stdev = sd(relab, na.rm = T)) -> relab
  
  out <- data.frame(cbind(prev, round(relab[,-1],1)), feature=mag, 
                    species=full.taxonomy$species[full.taxonomy$original_bin==mag])
  return(out)
}))

Bifid_summary2[Bifid_summary2$species != "Bifidobacterium_dentium",]
write.table(Bifid_summary2, "output_files/revised_MAG/Figure2/Bifido_summary.tsv")

Bifid_summary2[Bifid_summary2$species == "Bifidobacterium_longum",] -> blong
blong <- blong[order(blong$average),]
blong$relab <- blong$average/10^6

#############################
# Bacteroides vs Prevotella #
#############################
x <- full.taxonomy$species[full.taxonomy$feature%in%select]
prev <- x[grep("Prevotella",x)] # list in top 60
prev[prev %in% unique(data$species[(data$coef < 0 & data$metadata == "lifestyle" & data$qval<0.05)])] # assoc. w/ rural
prev[prev %in% unique(data$species[(data$coef > 0 & data$metadata == "lifestyle" & data$qval<0.05)])] # assoc. w/ urban

bact <- x[grep("Bacteroides|Phocaeicola",x)] # list in top 60
bact[bact %in% unique(data$species[(data$coef < 0 & data$metadata == "lifestyle" & data$qval<0.05)])] # assoc. w/ rural
bact[bact %in% unique(data$species[(data$coef > 0 & data$metadata == "lifestyle" & data$qval<0.05)])] # assoc. w/ urban


n <- c(full.taxonomy$original_bin[grep("Bacteroides|Phocaeicola",full.taxonomy$genus)],
       full.taxonomy$original_bin[grep("Prevotella",full.taxonomy$genus)])
corr.mat <- data.frame(input.df[,colnames(input.df) %in% n])
colnames(corr.mat) <- unlist(lapply(colnames(corr.mat), function(x) {
  x <- full.taxonomy$species[which(full.taxonomy$original_bin==x)]
  return(x)}))
psych::corr.test(corr.mat,method = "kendall") -> corr.out
rho <- corr.out$r
pval <- matrix(p.adjust(corr.out$p, method = "fdr"), ncol = ncol(corr.out$p), dimnames = dimnames(corr.out$p))
rho <- rho[rownames(rho) %in% x,colnames(rho) %in% x]
pval <- pval[rownames(pval) %in% x,colnames(pval) %in% x]
sel <- grep("Prevotella", rownames(rho))
sel <- grep("Prevotella", rownames(rho))
x <- round(rho[sel,-sel],2) # Rho
y <- round(pval[sel,-sel],4) # pval
matrix(paste(round(rho[sel,-sel],2), "|", round(pval[sel,-sel],4), sep = ""), ncol=3, dimnames = dimnames(x)) # top 60


################
## plot top 60 #
################
heat.data <- data[order(data$coef, decreasing = T),]
heat.data <- heat.data[heat.data$metadata=="lifestyle",]
heat.data <- heat.data[heat.data$slope!=0,]
heat.data <- heat.data[heat.data$qval<0.05,]
heat.data$abs.coef <- abs(heat.data$coef)

meanRel <- data.frame(original_bin=colnames(input.df[,-1]), meanRel=colMeans(input.df[,-1]/10^6), row.names = NULL)
meanRel <- merge.data.frame(meanRel, data.frame(feature=full.taxonomy$feature, original_bin=full.taxonomy$original_bin))
meanRel$original_bin <- NULL

heat.data <- merge.data.frame(heat.data,meanRel, by="feature")
heat.data$norm.coef <- (heat.data$abs.coef/heat.data$meanRel)*heat.data$slope
scale.base <- max(abs(heat.data$norm.coef))/3
heat.data$norm.coef <- heat.data$norm.coef/scale.base
heat.data$abs.norm.coef <- abs(heat.data$norm.coef)
heat.data <- heat.data[order(heat.data$abs.norm.coef),]
old.heat.data <- heat.data # save for later

summary(heat.data$coef)
summary(heat.data$norm.coef)

heat.data <- heat.data[heat.data$feature %in% select,]

phylcols <- read.delim("input_files/phylum_colours_bar.tsv")
phylum.cols <- phylcols$cols
names(phylum.cols) <- phylcols$phylum

hist(heat.data$coef)
hist(heat.data$norm.coef)

# order by coefficient
heat.data <- heat.data[order(heat.data$coef),]
rownames(heat.data) <- NULL
table(duplicated(heat.data$feature)) # merging check. all should be FALSE
heat.data$species

#mark taxon in bold
df.indomee <- df.indomee[df.indomee$qval<0.05,]
df.indomee <- df.indomee[df.indomee$value=="lifestyle",]
indomee_signif <- heat.data$feature[heat.data$feature %in% df.indomee$feature]
heat.data$xlabel <- as.character(heat.data$species)
heat.data$xlabel[heat.data$feature %in% indomee_signif] <- paste("**",heat.data$xlabel[heat.data$feature %in% indomee_signif],"**",sep = "") 

indomee_present <- input.df[!input.df$sampleid%in%metadata$sampleid[grep("SIN", metadata$population)],]
table(colSums(indomee_present[,-1]) > 0)

# order labels
xlabel.order <- as.character(heat.data$species)
indomee.species <- as.character(heat.data$species[heat.data$feature %in% indomee_signif])
xlabel.order[xlabel.order %in% indomee.species] <- paste("**",xlabel.order[xlabel.order %in% indomee.species],"**", sep="")
heat.data$xlabel <- factor(heat.data$xlabel, xlabel.order)

a <- ggplot(heat.data) + theme_classic() +
  scale_y_discrete(expand = c(0,0))+
  geom_tile(aes(y="phylum",x=xlabel, fill=phylum)) +
  scale_fill_manual(values=phylum.cols)+
  theme(plot.margin=margin(0,0,0,0,"cm"),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        legend.location = "plot")+
  ylab("Phylum")
#ggtitle("Top 50 features with significant associations (by coefficient)")

legend_a <- as_ggplot(cowplot::get_legend(a)) + theme(plot.margin=margin(0,0,0,0,"cm"))
panel_a <- as_ggplot(cowplot::get_panel(a)) + theme(plot.margin=margin(0,0,0,0,"cm"))

b <- ggplot(heat.data) + theme_classic() +
  scale_y_discrete(expand = c(0,0))+
  geom_tile(aes(y="A",x=xlabel, fill=coef))+
  scale_fill_gradient2(high="royalblue3", mid = "grey95", low="salmon",midpoint = 0,
                       breaks = c(-3,-2,-1, 0, 1,2,3))+
  theme(plot.margin=margin(0,0,0,0,"cm"),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        legend.location = "plot")+
  ylab("Raw Coefficient")
#ggtitle("Top 50 features with significant associations (by coefficient)")


legend_b <- as_ggplot(cowplot::get_legend(b)) + theme(plot.margin=margin(0,0,0,0,"cm"))
panel_b <- as_ggplot(cowplot::get_panel(b)) + theme(plot.margin=margin(0,0,0,0,"cm"))

c <- ggplot(heat.data) + theme_classic() +
  scale_y_discrete(expand = c(0,0))+
  geom_tile(aes(y="A",x=xlabel, fill=norm.coef))+
  scale_fill_gradient2(high="royalblue3", mid = "grey95", low="salmon",midpoint = 0,
                       breaks = c(-3,-2,-1, 0, 1,2,3))+
  theme(plot.margin=margin(0,0,0,0,"cm"),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        legend.location = "plot")+
  ylab("Normalised Coef.")


legend_c <- as_ggplot(cowplot::get_legend(c)) + theme(plot.margin=margin(0,0,0,0,"cm"))
panel_c <- as_ggplot(cowplot::get_panel(c)) + theme(plot.margin=margin(0,0,0,0,"cm"))

d <- ggplot(heat.data,aes(y="A",x=xlabel)) + theme_classic() +
  scale_y_discrete(expand = c(0,0)) +
  theme(plot.margin=margin(0,0,0,0,"cm"),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.x = element_markdown(angle=90, hjust=1, vjust=0.5),
        #axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        legend.location = "plot")+
  ylab("")

x.axis <- as_ggplot(cowplot::get_x_axis(d))

# draw plot
lay <- rbind(c(1,1,1,1,2),
             c(3,3,3,3,4),
             c(5,5,5,5,6),
             c(7,7,7,7,6),
             c(7,7,7,7,6))

grid.arrange(b+theme(legend.position = "none"),legend_b,
             c+theme(legend.position = "none"),legend_c,
             a+theme(legend.position = "none"),legend_a,
             x.axis+ theme(plot.margin=margin(0,0,0,.5,"cm")),
             layout_matrix=lay)

# save plot
ggsave(plot=grid.arrange(b+theme(legend.position = "none"),legend_b,
                         c+theme(legend.position = "none"),legend_c,
                         a+theme(legend.position = "none"),legend_a,
                         x.axis+ theme(plot.margin=margin(0,0,0,.5,"cm")),
                         layout_matrix=lay),
       filename = "output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/cont_p20/figures/Top60.svg",
       width = 55,height = 15, units = "cm",dpi = "print", device = "svg")


# Figure 2E
lay <- rbind(c(1,1,1,1,2),
             c(3,3,3,3,4),
             c(5,5,5,5,6),
             c(5,5,5,5,6),
             c(5,5,5,5,6))

grid.arrange(b+theme(legend.position = "none"),legend_b,
             c+theme(legend.position = "none"),legend_c,
             x.axis+ theme(plot.margin=margin(0,0,0,.5,"cm")),
             layout_matrix=lay)

lay <- rbind(c(1),c(2),c(3),c(3),c(3))
dir.create("figout/Figure2/",showWarnings = F)
ggsave(plot=grid.arrange(b+theme(legend.position = "none"),
             c+theme(legend.position = "none"),
             x.axis+ theme(plot.margin=margin(0,0,0,.5,"cm")),
             layout_matrix=lay),
       file="figout/Figure2/Fig2E_MAASLIN_SPMP.png", width = 10,height = 5,device = "png",dpi = 1200, scale = 1.2)

ggsave(plot=grid.arrange(b+theme(legend.position = "none"),
                         c+theme(legend.position = "none"),
                         x.axis+ theme(plot.margin=margin(0,0,0,.5,"cm")),
                         layout_matrix=lay),
       file="figout/Figure2/Fig2E_MAASLIN_SPMP.svg", width = 10,height = 5,device = "svg",dpi = 1200, scale = 1.2)

saveRDS(object = grid.arrange(b+theme(legend.position = "none"),
                              c+theme(legend.position = "none"),
                              x.axis+ theme(plot.margin=margin(0,0,0,.5,"cm")),
                              layout_matrix=lay),
        file = "figout/Figure2/Fig2E_MAASLIN_SPMP.rds")


###########
# Heatmap #
###########
library(reshape2)
metadata <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/input.metadata_p20.tsv", header=T)
pop.ord <- readRDS("input_files/pop.ord.rds")

data <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/cont_p20/all_results_MGYG_wTax.tsv", header=T)
data <- data[data$metadata=="lifestyle" &data$qval<0.05,]
signif.features <- unique(data$feature)

RA <- read.delim("input_files/MAGs/revised/Hybrid_wSPMP/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv")[,c("original_bin", unique(metadata$sampleid))]
RA <- merge.data.frame(full.taxonomy[,c("feature","genome","original_bin")], RA)
rownames(RA) <- RA$feature
RA <- RA[,unique(metadata$sampleid)]


RA <- t(RA)
RA <- data.frame(RA, sampleid=rownames(RA))
RA <- melt(RA, id.vars = "sampleid",variable.name = "feature",value.name = "RA")
RA <- merge.data.frame(RA, full.taxonomy[,c("feature","original_bin","phylum","class","family","genus","species")], by="feature")  
RA <- merge.data.frame(RA, metadata[,c("sampleid","population","lifestyle")], by="sampleid")
RA <- RA[RA$feature %in% signif.features,]
RA$population <- factor(RA$population, pop.ord)

# get coef and sign from maaslin
subdata <- data
subdata <- subdata[order(subdata$coef),]
coef.ord <- as.character(unique(subdata$feature))
RA$feature <- factor(RA$feature, levels=coef.ord)
RA <- RA[order(RA$feature),]
RA$species <- factor(RA$species, levels = unique(RA$species))

#################
# Selected Taxa #
#################
Bifi <- grep("Bifidobacterium",x = RA$species)
Trep <- grep("Treponema",x = RA$species)
Pse <- grep("Pseudomonadota",x = RA$phylum)
Prev <- grep("Prevotella",x = RA$genus)
Bact <- grep("Bacteroides|Phocaei",x = RA$genus)
Faecali <- grep("Faecalibacterium",x = RA$genus)
Blautia <- grep("Blautia",x = RA$genus)

heat.df <- RA[c(Bifi,Trep,Prev,Bact,Faecali,Blautia),]
heat.df[grep("Bacteroides|Phocaei", heat.df$genus),"genus"] <- "Bacteroides|Phocaeicola"
heat.df$genus <- factor(heat.df$genus, levels = c("Bifidobacterium", "Treponema_D","Prevotella","Bacteroides|Phocaeicola","Faecalibacterium", "Escherichia","Blautia"))
mid <- median(heat.df$RA)

# rearrange factors
heat.df <- heat.df[order(heat.df$feature),]
heat.df$species <- factor(heat.df$species, levels = unique(as.character(heat.df$species)))

p <- ggplot(heat.df, aes(x=sampleid,y=species))+ theme_minimal(base_size = 15)+
  geom_tile(aes(fill=-log10(RA)),width=1.05) +
  facet_grid(genus~population, space = "free", scales = "free")+
  scale_fill_gradient2(high = "grey20", mid="purple3", low = "darkorchid1",na.value = "grey20",midpoint = 2) +
  theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks.x = element_blank()) +
  theme(legend.position = "bottom")+
  theme(panel.background = element_blank(), panel.border = element_rect(linewidth = 0.5,colour = "white",fill = NA), panel.grid = element_blank())+
  theme(strip.text.y = element_blank(), strip.text.x = element_text(angle=90, hjust = 0, vjust = 0.5)) +
  theme(axis.text.y = element_text(size=8, hjust=1, vjust=0))+
  theme(panel.spacing.x = unit(-1, "mm"), panel.spacing.y = unit(-1, "mm"))

p

ggsave(p,
       filename = "output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/cont_p20/figures/RA_heatmap_p20_selected.svg",
       width = 55,height = 21, units = "cm",dpi = "print")

###############
# TOP Results #
###############
data$abs.coef <- abs(data$coef)
select <- select_taxon_by_coef(data,coef.column = "abs.coef",n = 100)

heat.df <- RA
heat.df <- RA[RA$feature%in%select,]
heat.df[grep("Bacteroides|Phocaei", heat.df$genus),"genus"] <- "Bacteroides|Phocaeicola"
mid <- median(heat.df$RA)

# get coef and sign from maaslin to order MAG
subdata <- data[data$feature%in%unique(heat.df$feature),]
subdata <- subdata[order(subdata$coef),]
coef.ord <- rev(as.character(subdata$feature))
heat.df$feature <- factor(heat.df$feature, levels = coef.ord)
heat.df <- heat.df[order(heat.df$feature),]
heat.df$species <- factor(heat.df$species, levels = unique(as.character(heat.df$species)))

# get average coef and sign from maaslin to order Genus
library(dplyr)
subdata[grep("Bacteroides|Phocaei", subdata$genus),"genus"] <- "Bacteroides|Phocaeicola"
subdata %>% group_by(genus) %>% summarise(sum.coef=sum(coef)) -> sumcoef.genus
sumcoef.genus <- sumcoef.genus[order(sumcoef.genus$sum.coef, decreasing = T),]
genus.ord <- sumcoef.genus$genus
heat.df$genus <- factor(heat.df$genus, genus.ord)

heat.df$sampleid <- factor(heat.df$sampleid, levels=metadata$sampleid)

p  <- ggplot(heat.df, aes(x=sampleid,y=species))+theme_minimal(base_size = 11)+
  geom_tile(aes(fill=-log10(RA)),width=1.05) +
  facet_grid(genus~population, space = "free", scales = "free")+
  scale_fill_gradient2(high = "grey20", mid="purple3", low = "darkorchid1",na.value = "grey20",midpoint = 2) +
  theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks.x = element_blank()) +
  theme(axis.text.y = element_blank())+
  theme(legend.position = "bottom")+
  #theme(panel.background = element_blank(), panel.border = element_rect(linewidth = 0.3,colour = "white",fill = NA), panel.grid = element_blank())+
  #theme(panel.spacing.x = unit(-1, "mm"), panel.spacing.y = unit(-1, "mm")) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid = element_blank())+
  theme(panel.spacing.x = unit(0.3, "mm"), panel.spacing.y = unit(0.3, "mm")) +
  theme(strip.text.x = element_text(angle=90, hjust = 0, vjust = 0.5), strip.text.y = element_text(angle=0, hjust = 0, vjust = 0.5)) +
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0))

p

ggsave(p,
       filename = "output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/cont_p20/figures/RA_heatmap_p20_selected_top100.svg",
       width = 55,height = 30, units = "cm",dpi = "print")
saveRDS(p, file="figout/Figure2/FigS2C_RA_heatmap_p20_top100.rds")

