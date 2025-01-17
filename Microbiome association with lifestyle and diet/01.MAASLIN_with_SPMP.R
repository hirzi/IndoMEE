library(ggplot2)
library(ggpubr)
library(gridExtra)

## custom function
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

# read input
input.df <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/input.comm_p20.tsv", header=T)
metadata <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/input.metadata_p20.tsv", header=T)
dim(input.df)
df <- merge.data.frame( metadata, input.df, by="sampleid")
rownames(df) <- df$sampleid

# load taxonomy
full.taxonomy <- read.csv("input_files/MAGs/revised/spore_taxonomy.csv", header=T)
rownames(full.taxonomy) <- full.taxonomy$original_bin
full.taxonomy <- full.taxonomy[full.taxonomy$original_bin%in%colnames(df),]
row.names(full.taxonomy) <- NULL
full.taxonomy$feature <- full.taxonomy$original_bin

# load results and filter for significance
data <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/cont_p20/lifestyle_results_wtax.tsv", header=T)
data$slope <- sign(data$coef)

signif.df <- data[data$qval<0.05,]
signif.df <- signif.df[order(abs(signif.df$coef), decreasing = T),]

length(unique(data$feature))
length(unique(signif.df$feature))
100*length(unique(signif.df$feature))/length(unique(data$feature))

# chisqr by signs
data[data$qval >0.05,"slope"] <- 0

table(data$phylum, data$slope)
chisq.test(table(data$phylum, data$slope))
chisq.test(table(data$family, data$slope))
chisq.test(table(data$genus, data$slope))

# bifido/treponema corr
n <- c(full.taxonomy$original_bin[grep("Bifidobacterium",full.taxonomy$genus)],
       full.taxonomy$original_bin[grep("Treponema",full.taxonomy$genus)])
corr.mat <- data.frame(input.df[,colnames(input.df) %in% n])
colnames(corr.mat) <- full.taxonomy$species[full.taxonomy$original_bin %in% n]
corr.mat <- corr.mat[,sort(colnames(corr.mat))]
psych::corr.test(corr.mat,method = "kendall") -> corr.out
corr.out$r[c(1:5),-c(1:5)]

# Bacteroides/Prevotella corr
n <- c(full.taxonomy$original_bin[grep("Bacteroides|Phocaei",full.taxonomy$genus)],
       full.taxonomy$original_bin[grep("Prevotella",full.taxonomy$genus)])
corr.mat <- data.frame(input.df[,colnames(input.df) %in% n])
colnames(corr.mat) <- full.taxonomy$species[full.taxonomy$original_bin %in% n]
corr.mat <- corr.mat[,sort(colnames(corr.mat))]
psych::corr.test(corr.mat,method = "kendall") -> corr.out
sel <- grep("Prevotella",colnames(corr.out$r))
corr.out$r[sel,-sel]

# volcano
data$enrichment <- ifelse(sign(data$coef)<0, "Rural Indonesia", "Singapore")
data[data$qval>0.05,"enrichment"] <- "not significant"

ggplot(data) + theme_bw() +
  geom_point(aes(x=coef, y=-log10(qval), col=enrichment))+
  scale_color_manual(values = c("grey","salmon","royalblue3"))+
  geom_hline(yintercept = -log10(0.05), col="red", size=1, linetype=2)

# top 30 heatmap
heat.data <- data[order(abs(data$coef), decreasing = T),]
heat.data <- heat.data[heat.data$slope!=0,]
select <- heat.data[1:100,"feature"]
heat.data <- heat.data[heat.data$feature %in% select,]
heat.data$x <- paste0(heat.data$species, "|", heat.data$feature)
heat.data$transition <- ifelse(heat.data$slope<0,"Rural Indonesian", "Singapore and Denpasar")

ggplot(heat.data) + theme_classic() +
  geom_tile(aes(y=transition,x=x, fill=coef))+
  scale_fill_gradient2(high="royalblue3", mid = "grey95", low="salmon")+
  theme(axis.text.x = element_text(angle=60, hjust=1, vjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(1,1,1,1, unit = "cm"))+
  ggtitle("Top 50 features with significant associations (by coefficient)")


####
#### Top 60
select_taxon_by_coef <- function(maaslin.output, coef.column, n){
  maaslin.output <- maaslin.output[order(maaslin.output[, coef.column], decreasing = T),]
  select <- maaslin.output$feature[1:n]
  select}

heat.data <- data[order(abs(data$coef), decreasing = T),]
heat.data <- heat.data[heat.data$slope!=0,]
heat.data <- heat.data[heat.data$qval<0.05,]
heat.data$abs.coef <- abs(heat.data$coef)

meanRel <- data.frame(feature=colnames(input.df[,-1]), meanRel=colMeans(input.df[,-1]/10^6), row.names = NULL)
heat.data <- merge.data.frame(heat.data,meanRel, by="feature")
heat.data$norm.coef <- (heat.data$abs.coef/heat.data$meanRel)*heat.data$slope
scale.base <- max(abs(heat.data$norm.coef))/3
heat.data$norm.coef <- heat.data$norm.coef/scale.base
heat.data$abs.norm.coef <- abs(heat.data$norm.coef)
heat.data <- heat.data[order(heat.data$abs.norm.coef),]
old.heat.data <- heat.data # save for later

summary(heat.data$coef)
summary(heat.data$norm.coef)

select <- select_taxon_by_coef(heat.data,coef.column = "abs.coef",n = 60)
saveRDS(select, file="output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/cont_p20/top60_features.rds")

heat.data <- heat.data[heat.data$feature %in% select,]
heat.data$transition <- ifelse(heat.data$slope<0,"Remote and Rural Indonesia", "Singapore and Denpasar")
table(heat.data$transition)

coef.order <- order(heat.data$coef)
heat.data$species <- factor(heat.data$species, heat.data$species[coef.order])

phylcols <- read.delim("input_files/phylum_colours_bar.tsv")
phylum.cols <- phylcols$cols
names(phylum.cols) <- phylcols$phylum

hist(heat.data$coef)
hist(heat.data$norm.coef)

a <- ggplot(heat.data) + theme_classic() +
  scale_y_discrete(expand = c(0,0))+
  geom_tile(aes(y="phylum",x=species, fill=phylum)) +
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
  geom_tile(aes(y="A",x=species, fill=coef))+
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
  geom_tile(aes(y="A",x=species, fill=norm.coef))+
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

d <- ggplot(heat.data,aes(y="A",x=species)) + theme_classic() +
  scale_y_discrete(expand = c(0,0)) +
  theme(plot.margin=margin(0,0,0,0,"cm"),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
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
grid.arrange(b+theme(legend.position = "none"),
             c+theme(legend.position = "none"),
             x.axis+ theme(plot.margin=margin(0,0,0,.5,"cm")),
             layout_matrix=lay)

saveRDS(object = grid.arrange(b+theme(legend.position = "none"),
                              c+theme(legend.position = "none"),
                              x.axis+ theme(plot.margin=margin(0,0,0,.5,"cm")),
                              layout_matrix=lay),
        file = "output_files/revised_MAG/Figure2/Fig2E_MAASLIN_SPMP.rds")


###########
# Heatmap #
###########
heat.df <- RA[c(Bifi,Trep,Prev,Bact,Faecali,Blautia),]
heat.df[grep("Bacteroides|Phocaei", heat.df$genus),"genus"] <- "Bacteroides|Phocaeicola"
heat.df$genus <- factor(heat.df$genus, levels = c("Bifidobacterium", "Treponema_D","Prevotella","Bacteroides|Phocaeicola","Faecalibacterium", "Escherichia","Blautia"))
mid <- median(heat.df$RA)

# get coef and sign from maaslin
subdata <- data[data$feature%in%unique(heat.df$original_bin),]
subdata <- subdata[order(subdata$coef),]
coef.ord <- rev(as.character(subdata$feature))
heat.df$original_bin <- factor(heat.df$original_bin, levels = coef.ord)
heat.df <- heat.df[order(heat.df$original_bin),]
heat.df$species <- factor(heat.df$species, levels = unique(as.character(heat.df$species)))

ggplot(heat.df, aes(x=sampleid,y=species))+theme_minimal(base_size = 15)+
  geom_tile(aes(fill=-log10(RA)),width=1.05) +
  facet_grid(genus~population, space = "free", scales = "free")+
  scale_fill_gradient2(high = "grey20", mid="purple3", low = "darkorchid1",na.value = "grey20",midpoint = 2) +
  theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks.x = element_blank()) +
  theme(legend.position = "bottom")+
  theme(panel.background = element_blank(), panel.border = element_rect(linewidth = 0.5,colour = "white",fill = NA), panel.grid = element_blank())+
  theme(strip.text.y = element_blank(), strip.text.x = element_text(angle=90, hjust = 0, vjust = 0.5)) +
  theme(panel.spacing.x = unit(-1, "mm"), panel.spacing.y = unit(-1, "mm")) -> p

p

ggsave(p,
        filename = "output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/cont_p20/figures/RA_heatmap_p20.svg",
       width = 55,height = 21, units = "cm",dpi = "print")

# TOP Results
select <- select_taxon_by_coef(old.heat.data,coef.column = "abs.coef",n = 100)

heat.df <- RA
heat.df <- RA[RA$original_bin%in%select,]
heat.df[grep("Bacteroides|Phocaei", heat.df$genus),"genus"] <- "Bacteroides|Phocaeicola"
mid <- median(heat.df$RA)

# get coef and sign from maaslin to order MAG
subdata <- data[data$feature%in%unique(heat.df$original_bin),]
subdata <- subdata[order(subdata$coef),]
coef.ord <- rev(as.character(subdata$feature))
heat.df$original_bin <- factor(heat.df$original_bin, levels = coef.ord)
heat.df <- heat.df[order(heat.df$original_bin),]
heat.df$species <- factor(heat.df$species, levels = unique(as.character(heat.df$species)))

# get average coef and sign from maaslin to order Genus
library(dplyr)
subdata[grep("Bacteroides|Phocaei", subdata$genus),"genus"] <- "Bacteroides|Phocaeicola"
subdata %>% group_by(genus) %>% summarise(sum.coef=sum(coef)) -> sumcoef.genus
sumcoef.genus <- sumcoef.genus[order(sumcoef.genus$sum.coef, decreasing = T),]
genus.ord <- sumcoef.genus$genus
heat.df$genus <- factor(heat.df$genus, genus.ord)

heat.df$sampleid <- factor(heat.df$sampleid, levels=metadata$sampleid)

ggplot(heat.df, aes(x=sampleid,y=species))+theme_minimal(base_size = 11)+
  geom_tile(aes(fill=-log10(RA)),width=1.05) +
  facet_grid(genus~population, space = "free", scales = "free")+
  scale_fill_gradient2(high = "grey20", mid="purple3", low = "darkorchid1",na.value = "grey20",midpoint = 2) +
  theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks.x = element_blank()) +
  theme(axis.text.y = element_blank())+
  theme(legend.position = "bottom")+
  #theme(panel.background = element_blank(), panel.border = element_rect(linewidth = 0.3,colour = "white",fill = NA), panel.grid = element_blank())+
  #theme(panel.spacing.x = unit(-1, "mm"), panel.spacing.y = unit(-1, "mm")) +
  theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid = element_blank())+
  theme(panel.spacing.x = unit(0.25, "mm"), panel.spacing.y = unit(0.25, "mm")) +
  theme(strip.text.x = element_text(angle=90, hjust = 0, vjust = 0.5), strip.text.y = element_text(angle=0, hjust = 0, vjust = 0.5)) +
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) -> p

p

ggsave(p,
       filename = "output_files/revised_MAG/MAASLIN2/maaslin2_SPMP_CLR_LsGrd/cont_p20/figures/RA_heatmap_p20_top100.svg",
       width = 55,height = 30, units = "cm",dpi = "print")
saveRDS(p, file="output_files/revised_MAG/FigureS2/FigS2B_RA_heatmap_p20.rds")
