# setup
library(ggplot2)
library(gridExtra)
library(vegan)

# custom functions
convert_to_CPM <- function(comm.df){
  comm.df <- comm.df*10^6
  comm.df}
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
axis_labels <- function(cap){
  eig <- cap$CA$eig
  var <- eig/sum(eig)
  axis_labs <- paste("PCo",1:length(names(var)), " (", round(var*100,1),"%)", sep="")
  axis_labs
}
lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  l <- list(a = format(unname(coef(m)[1]), digits = 2),
            b = format(unname(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3))
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  as.character(as.expression(eq))
}


# load popid
mypop <- readRDS("output_files/clindata_imputed.rds")[,c("sampleid","population")]
sampleid <- unique(as.character(mypop$sampleid))
removed_samples <- "B1|B2|B3|REP|5B|DG|FC|POS|NTC|EB|Sum|X|JKT|RPN001"
sampleid <- sampleid[-grep(removed_samples, sampleid)]

pop.ord <- c("Asmat_DAI","Punan_PBS","Punan_HUL","Punan_SUL","Punan_RPN","Basap_BRU","Balinese_PDW","Balinese_DPS")
mypop <- mypop[mypop$population %in% pop.ord,]
mypop$population <- factor(mypop$population, pop.ord)
mypop <- mypop[order(mypop$population),]
sample.ord <- as.character(unique(mypop$sampleid))

# load popcols
mycols <- read.delim("input_files/MAGs/MAGs_pop_mycols.tsv", header=T)
mycols <- subset(mycols, group %in% pop.ord)
tmp <- mycols$mycols ; names(tmp) <- mycols$group
mycols <- tmp

# load diet data
dat <- read.table("input_files/final_metadata_ClinLs.tsv")
select.var <- c("Rice","Sago","Pork","Chicken","Fish","Egg","Milk","Hunting","Watching_TV")
pcdata <- dat[,c("sampleid","population",select.var)]
pcdata$population <- factor(pcdata$population, pop.ord)
pcdata$sampleid <- factor(pcdata$sampleid, sample.ord)
pcdata <- pcdata[order(pcdata$sampleid),]
rownames(pcdata) <- pcdata$sampleid
pcdata <- na.omit(pcdata)
pcdata <- pcdata[,-c(1,2)]
diet.scaled <- apply(pcdata, 2, function(x) (x-mean(x))/sd(x))

# old pca data
# diet.scaled2 <- data.frame(readRDS("output_files/MAGs/pmv.diet.scaled.rds"))
# diet.scaled2 <- diet.scaled2[rownames(diet.scaled),]
# table(diet.scaled2 == diet.scaled)

# diet pca
diet.pca <- rda(diet.scaled)
diet.axes <- gsub(pattern = "Co",replacement = "C", axis_labels(diet.pca))

myaxes <- c(1,2)
summary.pca <- summary(diet.pca)


#ggdata <- data.frame(sampleid=rownames(diet.pca$CA$u), diet.pca$CA$u[,myaxes],row.names = NULL)
# ggdata <- data.frame(sampleid=rownames(summary.pca$sites), summary.pca$sites[,myaxes],row.names = NULL) # old version
ggdata <- data.frame(sampleid=rownames(scores(diet.pca,choices = myaxes)$sites), scores(diet.pca,choices = myaxes)$sites,row.names = NULL)
ggdata <- merge.data.frame(ggdata,mypop, by="sampleid")
ggdata$population <- factor(ggdata$population,levels = pop.ord)
ggdata$sampleid <- factor(ggdata$sampleid, sample.ord)
ggdata <- ggdata[order(ggdata$sampleid),]
ggdata$urbanisation <- ifelse(ggdata$population %in% c("Asmat_DAI", "Punan_HUL"), "remote",
                              ifelse(ggdata$population %in% c("Basap_BRU", "Balinese_PDW"), "rural","urban"))

xtitle <- diet.axes[myaxes[1]]; ytitle <- diet.axes[myaxes[2]]
breaks <- seq(-2,2,by=0.1)

dim(ggdata)

baseplot <- ggplot() +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
        legend.position = "bottom",
        legend.title=element_blank(),
        legend.text = element_text(size=12),
        legend.key.size = unit(6,"mm"),
        axis.text=element_blank(),
        axis.ticks = element_blank())

 # fit env
diet.fit <- envfit(diet.pca, diet.scaled, permutations=9999)
print(diet.fit)

diet.vec <- diet.fit$vectors
diet.vec <- data.frame(diet.vec$arrows)/1
pvalselect <- diet.fit$vectors$pvals <1
diet.vec <- diet.vec[pvalselect,]

limit_x <- c(min(-1*ggdata$PC1)*1.1, max(-1*ggdata$PC1)*2) ; limit_y <- c(min(-1*ggdata$PC2)*2.5, max(-1*ggdata$PC2)*1.1)

a <- baseplot + 
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept=0, lty="dashed") +
  geom_point(data=ggdata, aes(x=-1*PC1, y=-1*PC2, fill=population, 
                              shape=urbanisation, size=urbanisation),
             col="black") +
  scale_shape_manual(values = c(21,22,24))+
  scale_size_manual(values = c(5,4.7,4.5))+
  scale_fill_manual(values=mycols) +
  scale_x_continuous(expand = c(0,0),
                     breaks = breaks,
                     limits = limit_x)+
  scale_y_continuous(expand = c(0,0),
                     limits = limit_y,
                     breaks = breaks,
                     position = "left") +
  xlab(xtitle) + ylab(ytitle)+
  geom_segment(data = diet.vec, 
               aes(x = 0, y = 0, xend = -1*(PC1), yend = -1*(PC2)), 
               arrow = arrow(length = unit(1, "picas")),
               color = "red", size=1.5) +
  geom_text (data=diet.vec, aes(x=-1.3*PC1,y=-1.3*PC2,label=rownames(diet.vec)), 
             col="red", hjust=0.5, size=4,
             nudge_x = diet.vec$PC1/8,
             nudge_y = diet.vec$PC2/8) +
  labs(title = "PCA Lifestyle data", caption = paste0("Sample #:",nrow(diet.scaled)))

a

outdir="output_files/revised_MAG/revison_fig/Figure3/"
saveRDS(grid.arrange(a + theme(legend.position = "none")), file.path(outdir,"Fig3A.rds"))
ggsave(plot=a, filename = file.path(outdir,"Fig3A.png"),height = 5,width = 5,device = "png", units = "in",scale = 2)
ggsave(plot=a, filename = file.path(outdir,"Fig3A.svg"),height = 5,width = 5,device = "svg", units = "in",scale = 2)
ggsave(plot=a, filename = file.path(outdir,"Fig3A.pdf"),height = 5,width = 5,device = "svg", units = "in",scale = 2)


# Linear regression with microbes data
# Strong Filter
df <- read.delim("input_files/MAGs/revised/Hybrid/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv", header = T)
reads <- df[,colnames(df) %in% rownames(diet.scaled)]
# Ob <- readRDS("input_files/Ob_samples.rds")
# reads <- reads[,!colnames(reads) %in% Ob]
rownames(reads) <- df$original_bin
reads <- convert_to_CPM(reads)
reads <- prevalence_filter(reads, n=2)
reads <- as.matrix(t(reads))
dim(reads)

# select samples
comm.df <- reads[rownames(reads)%in%rownames(diet.scaled),]
comm.df <- comm.df[,colSums(comm.df)>1]
dim(comm.df)

# PCOA
comm.dist <- vegdist(comm.df+1, method="aitchison")
cap <- capscale(comm.dist~1)
pco.axes <- paste("Microbiome", axis_labels(cap))
pco_coord <- data.frame(sampleid=rownames(cap$CA$u), cap$CA$u)

# diet PCA
n65 <- rownames(comm.df)
# saveRDS(object = n65, "input_files/n65_sampleid.rds")
diet65 <- diet.scaled[rownames(diet.scaled) %in% n65,]
dim(diet65)
diet.pca <- rda(diet65)
diet.axes <- paste("Diet", gsub(pattern = "Co",replacement = "C", axis_labels(diet.pca)))
diet_coord <- data.frame(sampleid=rownames(diet.pca$CA$u), diet.pca$CA$u)

ggdata65 <- merge.data.frame(pco_coord[,c(1:3)], diet_coord[,c(1:3)], by="sampleid")
ggdata65 <- merge.data.frame(ggdata65, mypop, by="sampleid")
dim(ggdata65)

ggdata65$urbanisation <- ifelse(ggdata65$population %in% c("Asmat_DAI", "Punan_HUL"), "remote",
                                ifelse(ggdata65$population %in% c("Basap_BRU", "Balinese_PDW"), "rural","urban"))


b <- baseplot +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")+
  geom_point(data=ggdata65, aes(x=-1*PC1, y=MDS1, fill=population,
                                shape=urbanisation, size=urbanisation), col="black") +
  scale_shape_manual(values = c(21,22,24))+
  scale_size_manual(values = c(5,4.7,4.5))+
  scale_fill_manual(values=mycols) +
  xlab(diet.axes[1]) + ylab(pco.axes[1]) +
  geom_smooth(data=ggdata65, aes(x=-1*PC1, y=MDS1), method=lm) +
  annotate("text", x=(min(-1*ggdata65$PC1)+max(-1*ggdata65$PC1))/2,y=min(ggdata65$MDS1)*1.1, 
           label=lm_eqn(data.frame(x=-1*ggdata65$PC1, y=ggdata65$MDS1)),
           parse = TRUE, size=5)+
  labs(title = paste0("Linear Regression (n=", nrow(ggdata65),")"),
       caption = paste0("PCA:", paste0(colnames(diet65),collapse = "+")))

b

outdir="output_files/revised_MAG/revison_fig/Figure3/"
saveRDS(b+theme(legend.position = "none"), file.path(outdir,"Fig3B.rds"))
ggsave(plot=b, filename = file.path(outdir,"Fig3B.png"),height = 5,width = 5,device = "png", units = "in",scale = 2)
ggsave(plot=b, filename = file.path(outdir,"Fig3B.svg"),height = 5,width = 5,device = "svg", units = "in",scale = 2)
ggsave(plot=b, filename = file.path(outdir,"Fig3B.pdf"),height = 5,width = 5,device = "svg", units = "in",scale = 2)

# univariate regression#
library(lme4)
library(lmerTest)
fit <- lm(MDS1 ~ PC1, ggdata65)
summary(fit)

fitlme <- lmerTest::lmer(MDS1 ~ PC1 + (1|population), ggdata65)
summary(fitlme) -> out
out
anova(fitlme)
confint(fitlme)
fitlme

X <- model.matrix(fitlme)
beta <- fixef(fitlme)
eta <- as.matrix(X) %*% beta
var_fixed <- var(as.numeric(eta))

ranef.var <- data.frame(VarCorr(fitlme))[,c("grp","var1","vcov","sdcor")]
varex <- ranef.var$vcov
names(varex) <- ranef.var$grp
var_random <- varex["population"]
var_residual <- varex["Residual"]

# proportion of variance
var_total <- var_fixed + var_random + var_residual

100*var_random/var_total

########################################
# Between vs Within Population Compare #
########################################
sid <- mypop$sampleid
sid <- sid[-grep("MALRPN001",sid)]
Ob <- readRDS("input_files/Ob_samples.rds")


ls.df <- diet.scaled[rownames(diet.scaled)%in%sid,]
ls.scaled.pos <- apply(ls.df, 2, function(x) (x - min(x)) / (max(x) - min(x)))
lifestyle_chidist <- vegdist(ls.scaled.pos, method = "euclidean")
dist_df <- as.data.frame(as.table(as.matrix(lifestyle_chidist)),stringsAsFactors = F)
colnames(dist_df) <- c("Sample1","Sample2","Lifestyle_Dist")
dist_df <- subset(dist_df, dist_df$Sample1 < dist_df$Sample2)
dist_df <- dist_df[dist_df$Sample1!=dist_df$Sample2,]
dist_df$Pairs <- paste(dist_df$Sample1,dist_df$Sample2, sep="_")
plot_data <- dist_df

microb.comm.df <- data.frame(df[,colnames(df)%in%sid], row.names = df$original_bin)
microb.comm.df <- convert_to_CPM(microb.comm.df)
microb.comm.df <- prevalence_filter(microb.comm.df, n=2)
microb.comm.df <- as.matrix(t(microb.comm.df))
dim(microb.comm.df)
microb.dist <- vegdist(microb.comm.df+1, method = "aitchison")
dist_df <- as.data.frame(as.table(as.matrix(microb.dist)),stringsAsFactors = F)
colnames(dist_df) <- c("Sample1","Sample2","Microb_Dist")
dist_df <- subset(dist_df, dist_df$Sample1 < dist_df$Sample2)
dist_df <- dist_df[dist_df$Sample1!=dist_df$Sample2,]
dist_df$Pairs <- paste(dist_df$Sample1,dist_df$Sample2, sep="_")

countsample <- rownames(microb.comm.df)
countsample <- countsample[countsample %in% rownames(ls.scaled.pos)]
countsample <- gsub("BA0", "BALDPS0", countsample)
countsample <- gsub("BAPDW", "BALPDW", countsample)
countsample <- substr(countsample, 1,6)
table(countsample)

plot_data <- merge.data.frame(dist_df, plot_data[,c("Lifestyle_Dist","Pairs")], by="Pairs")

plot_data$Sample1 <- gsub(pattern = "BA0", "BALDPS0", plot_data$Sample1)
plot_data$Sample2 <- gsub(pattern = "BA0", "BALDPS0", plot_data$Sample2)
plot_data$Sample1 <- gsub(pattern = "BAPDW", "BALPDW", plot_data$Sample1)
plot_data$Sample2 <- gsub(pattern = "BAPDW", "BALPDW", plot_data$Sample2)

plot_data$cPairs <- ifelse(substr(plot_data$Sample1,start = 1,stop = 6) == substr(plot_data$Sample2,start = 1,stop = 6), "Within_Population", "Between_Population")
plot_data$Obese <- "no_obese"
plot_data$Obese [grep(pattern = paste0(Ob, collapse = "|"), plot_data$Pairs)] <- "with_obese"
plot_data$Obese [plot_data$Microb_Dist<65] <- "outlier"

ggplot()+ theme_classic(base_size=18)+
  geom_density(data=plot_data, aes(x=Microb_Dist, colour=cPairs), linewidth=1.5)+
  scale_x_continuous(breaks = seq(0,max(plot_data$Microb_Dist), by=10))+
  theme(legend.position = "top") +
  scale_colour_manual(values = c("royalblue3","brown"))
ggplot()+ theme_classic(base_size=18)+
  geom_density(data=plot_data, aes(x=Lifestyle_Dist, colour=cPairs),linewidth=1.5)+
  scale_x_continuous(breaks = seq(0,max(plot_data$Lifestyle_Dist), by=0.5))+
  theme(legend.position = "top") +
  scale_colour_manual(values = c("royalblue3","brown"))

outl <-subset(plot_data, Obese=="outlier")

c <-  baseplot +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")+
  geom_point(data=subset(plot_data, Obese!="outlier"), aes(x=Lifestyle_Dist, y=Microb_Dist, colour=cPairs, shape=Obese, alpha=cPairs), size=2)+
  scale_shape_manual(values = c(16,1))+
  scale_alpha_manual(values = c(0.35,0.30))+
  geom_smooth(data=subset(plot_data, Obese!="outlier"), aes(x=Lifestyle_Dist, y=Microb_Dist, col=cPairs), method = "lm", linewidth=1.8,se = F) +
  scale_colour_manual(values = c("royalblue","brown")) +
  ylab("Microbiome Distance (Aitchison)") + xlab("Lifestyle Distance (Euclidean)")+
  labs(title = paste0("Microbiome~Lifestyle association for Between vs Within Population"),
       caption = paste(length(unique(plot_data$Pairs))-nrow(outl),"pairwise relation across",length(unique(plot_data$Sample1)), "samples.", nrow(outl), "outlier points excluded."))

c

outdir="output_files/revised_MAG/revison_fig/FigureS3"
ggsave(plot=c, filename = file.path(outdir,"Fig3C.png"),height = 5,width = 4.8,device = "png", units = "in",scale = 2)
ggsave(plot=c, filename = file.path(outdir,"Fig3C.svg"),height = 5,width = 4.8,device = "svg", units = "in",scale = 2)
ggsave(plot=c, filename = file.path(outdir,"Fig3C.pdf"),height = 5,width = 4.8,device = "pdf", units = "in",scale = 2)
saveRDS(c, file.path(outdir,"Fig3C.rds"))

print(outl)
write.table(outl, "output_files/revised_MAG/revison_fig/additional_analyses/outlier_microbe_vs_ls_dist.tsv",row.names = F,quote = F,sep = "\t")

table(plot_data$cPairs)
length(plot_data$cPairs)

# Mantel test
x <- as.matrix(microb.dist)
y <- as.matrix(lifestyle_chidist)

# extend code
extend_sid_code <- function(sid){
  tmp <- sid
  n <- grep("BAPDW",tmp)
  tmp[n] <- gsub("BAPDW","BALPDW",tmp[n])
  n <- grep("BA0",tmp)
  tmp[n] <- gsub("BA0","BALDPS0",tmp[n])
  n <- grep("PBSHF",tmp)
  tmp[n] <- gsub("PBSHF","PBSHFE",tmp[n])
  return(tmp)
}

colnames(x) <- extend_sid_code(colnames(x))
rownames(x) <- extend_sid_code(rownames(x))
colnames(y) <- extend_sid_code(colnames(y))
rownames(y) <- extend_sid_code(rownames(y))

sid_match <- unique(colnames(x)[colnames(x) %in% colnames(y)], colnames(y)[colnames(y) %in% colnames(x)])
x <- as.dist(x[sid_match, sid_match])
y <- as.dist(y[sid_match, sid_match])

mask_group <- function(dist, group=c("between","within")){
  tmp <- as.matrix(dist)
  groups <- substr(colnames(tmp),1,6) # the first 6 character in a sampleid
  group_matrix <- outer(groups, groups, function(a, b) as.numeric(a != b))
  # between/within group matrix: 0 if same population (within), 1 if different population (between)
  
  if(group=="between"){
    tmp[group_matrix==0] <- NA
  } else if(group == "within"){
    tmp[group_matrix==1] <- NA
  }
  tmp <- as.dist(tmp)
  return(tmp)
}
x.btw <- mask_group(x,group = "between")
y.btw <- mask_group(y,group = "between")

x.wth <- mask_group(x, group = "within")
y.wth <- mask_group(y, group = "within")

mantel(x.btw, y.btw, permutations = 9999, na.rm=T)
mantel(x.wth, y.wth, permutations = 9999, na.rm=T)

################
# Mantel Tests #
################
library(ggpubr)
# population-level mentel
dirpath <- "input_files/population_level_distance/"
micro_dist <- read.csv(file.path(dirpath,"population_comm_distance_total.tsv"), sep = "\t", row.names = 1)
lifestyle_dist <- read.csv(file.path(dirpath,"lifestyle_comm_distance_average.tsv"), sep = "\t", row.names = 1)
geo_dist <- read.csv(file.path(dirpath,"Population_geoCoords_dist_metres.csv"), row.names = 1)

pop.ord <- colnames(micro_dist)[c(1,4,8,7,6,2,5,3,9,10,11)]
micro_dist <- micro_dist[pop.ord,pop.ord]
lifestyle_dist <- lifestyle_dist[pop.ord[pop.ord%in%colnames(lifestyle_dist)],pop.ord[pop.ord%in%colnames(lifestyle_dist)]]
geo_dist <- geo_dist[pop.ord, pop.ord]
indoMEE_pop <- pop.ord[pop.ord%in%colnames(lifestyle_dist)]

dist_df <- cbind(geo_dist = data.frame(as.dist(geo_dist[indoMEE_pop,indoMEE_pop]))/1000, 
                 micro_dist = data.frame(as.dist(micro_dist[indoMEE_pop,indoMEE_pop])), 
                 lifestyle_dist = data.frame(as.dist(lifestyle_dist[indoMEE_pop,indoMEE_pop])))
colnames(dist_df) <-  c("geo_dist","micro_dist","lifestyle_dist")

library(patchwork)

plot1 <- ggplot(dist_df, aes(x = geo_dist, y = micro_dist)) + 
  geom_point(size=5) + 
  geom_smooth(method = "lm") + 
  scale_y_continuous(limits = c(50, 150), breaks = seq(50, 150, 25)) +
  labs(x = "Distance (km)", y = "Microbiome divergence") + 
  theme_pubr(base_size = 18) + 
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major.y = element_line(linetype = 2))

plot2 <- ggplot(dist_df, aes(x = lifestyle_dist, y = micro_dist)) + 
  geom_point(size=5) + 
  geom_smooth(method = "lm") + 
  scale_y_continuous(limits = c(50, 150), breaks = seq(50, 150, 25)) +
  labs(x = "Lifestyle divergence", y = "Microbiome divergence") + 
  theme_pubr(base_size = 18) + 
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major.y = element_line(linetype = 2))

plot3 <- ggplot(dist_df, aes(x = geo_dist, y = lifestyle_dist)) + 
  geom_point(size=5) + 
  geom_smooth(method = "lm") + 
  #scale_y_continuous(limits = c(50, 150), breaks = seq(50, 150, 25)) +
  labs(x = "Distance (km)", y = "Lifestyle similarity") + 
  theme_pubr(base_size = 18) + 
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major.y = element_line(linetype = 2))

plot3 <- plot1 + plot2

plot3

# outdir="output_files/revised_MAG/FigureS2/"
outdir="output_files/revised_MAG/revison_fig/FigureS2"
dir.create(outdir)
saveRDS(plot3, file.path(outdir,"FigS2E.rds"))
ggsave(plot=plot3, filename = file.path(outdir,"FigS2E.png"),height = 5,width = 10,device = "png", units = "in",scale = 1)
ggsave(plot=plot3, filename = file.path(outdir,"FigS2E.svg"),height = 5,width = 10,device = "svg", units = "in",scale = 1)


################################
# Extract Factor Contributions #
################################
library(reshape2)
library(RColorBrewer)
library(ggpubr)

n65 <- rownames(comm.df)
# saveRDS(object = n65, "input_files/n65_sampleid.rds")
diet65 <- diet.scaled[rownames(diet.scaled) %in% n65,]
dim(diet65)
diet.pca <- rda(diet65)

good <- data.frame(goodness(diet.pca, choices = 1:9, display="species", model="CA"))

good2 <- t(apply(good, 1, function(x) diff(x)))
good2 <- data.frame(factor=rownames(good2),
                    PC1 = good$PC1,
                    good2)

good2.long <- melt(good2, variable.name = "PC",id.vars = "factor",value.name = "contribution")
factor.order <- rownames(good)[order(good$PC1)]
good2.long$factor <- factor(good2.long$factor, factor.order)

set.seed(9)
pc.cols=sample(c(brewer.pal(8,"Dark2"), brewer.pal(12,"Set3")), size=length(unique(good2.long$PC)))
p <- ggplot(good2.long) + theme_minimal(base_size = 16) +
  theme(legend.position = "bottom")+
  theme(panel.grid = element_blank())+
  geom_bar(aes(x=factor, y=100*contribution, fill=PC), stat="identity", position = "stack") +
  ylab("Factor Contribution (%)") +
  scale_fill_manual(values = pc.cols) +
  coord_flip() +
  labs(title="PCA - factor contribution", 
       subtitle = paste0("(n=", nrow(diet65), ", complete cases only)")) +
  guides(fill = guide_legend(title = "", nrow = 1))

p 

outdir="output_files/revised_MAG/revison_fig/FigureS3"
saveRDS(p, file.path(outdir,"FigS3C.rds"))
ggsave(plot=p, filename = file.path(outdir,"FigS3C.png"),height = 5,width = 10,device = "png", units = "in",scale = 0.8)
ggsave(plot=p, filename = file.path(outdir,"FigS3C.svg"),height = 5,width = 10,device = "svg", units = "in",scale = 0.8)
ggsave(plot=p, filename = file.path(outdir,"FigS3C.pdf"),height = 5,width = 10,device = "pdf", units = "in",scale = 0.8)


### By Population (complete cases only, n=73)
selectpop <- c("Punan_HUL","Basap_BRU","Balinese_PDW","Balinese_DPS")
# n65 <- readRDS("input_files/n65_sampleid.rds")
n65 <- rownames(comm.df)
mypop <- subset(mypop, sampleid %in% colnames(df))
mypop65 <- subset(mypop, sampleid %in% n65)
mypop65 <- subset(mypop, population %in% selectpop)
diet65 <- diet.scaled[rownames(diet.scaled) %in% n65,]
diet65 <- diet65[n65,]
comm.df <- comm.df[rownames(comm.df) %in% n65, ]
comm.df <- comm.df[n65,]
dim(comm.df)

comm.df.bypop <- list()
for(i in unique(mypop65$population)){
  tmp <- comm.df[rownames(comm.df) %in% mypop$sampleid[mypop$population==i],]
  tmp <- tmp[,colSums(tmp)>0]
  len <- length(comm.df.bypop)
  myname <- c(names(comm.df.bypop),i)
  comm.df.bypop[[len+1]] <- tmp
  names(comm.df.bypop)  <- myname
}

comm.dist.bypop.n65 <- list()
for(pop in names(comm.df.bypop)) {
  x <- comm.df.bypop[[pop]]
  x.dist <- vegdist(x+1, method="aitchison")
  len <- length(comm.dist.bypop.n65)
  nm <- c(names(comm.dist.bypop.n65),pop)
  comm.dist.bypop.n65[[len+1]] <- x.dist
  names(comm.dist.bypop.n65) <- nm
}

# PCoA
cap.bypop.n65 <- list()
cap.summary.bypop.n65 <- list()

for(pop in names(comm.dist.bypop.n65)){
  x <- as.dist(comm.dist.bypop.n65[[pop]])
  cap <- capscale(x~1)
  
  len <- length(cap.bypop.n65)
  nm <- c(names(cap.bypop.n65),pop)
  cap.bypop.n65[[len+1]] <- cap
  names(cap.bypop.n65) <- nm
  
  cap.summary <- summary(cap)
  cap.summary.bypop.n65[[len+1]] <- cap.summary
  names(cap.summary.bypop.n65) <- nm
  
}

# Diet PCA
diet.pca.bypop <- list()
diet.pca.summary.bypop <- list()
diet.envfit.bypop <- list()
dim(diet65)
for(pop in unique(mypop65$population)){
  sid <- as.character(mypop65$sampleid[mypop65$population==pop])
  x <- diet65[rownames(diet65)%in%sid,]
  dpca <- rda(x)
  
  len <- length(diet.pca.bypop)
  nm <- c(names(diet.pca.bypop),pop)
  diet.pca.bypop[[len+1]] <- dpca
  names(diet.pca.bypop) <- nm
  
  dpca.summary <- summary(dpca)
  diet.pca.summary.bypop[[len+1]] <- dpca.summary
  names(diet.pca.bypop) <- nm
  
  diet.fit <- envfit(ord=dpca, env=x, permutations=999)
  diet.envfit.bypop[[len+1]] <- diet.fit
  names(diet.envfit.bypop) <- nm
}


# draw lm plot
lm_plot_out <- list()
for(pop in unique(mypop65$population)){
  pco <- cap.bypop.n65[[pop]]
  ylab <- gsub("PCo", "Microbiome PCo", axis_labels(pco)[1])
  pco <- pco$CA$u
  pco <- data.frame(sampleid=rownames(pco), MDS1=pco[,1])
  
  dpca <- diet.pca.bypop[[pop]]
  xlab <- gsub("PCo","Diet PC", axis_labels(dpca)[1])
  dpca <- dpca$CA$u
  dpca <- data.frame(sampleid=rownames(dpca), PC1=dpca[,1])
  
  lmdata <- merge.data.frame(pco,dpca, by="sampleid")
  lmdata <- merge.data.frame(mypop65, lmdata, by="sampleid")
  
  
  out <- baseplot +
    geom_vline(xintercept = 0,linetype=2)+
    geom_hline(yintercept = 0,linetype=2)+
    theme(panel.grid = element_blank(), title = element_text(size = 14),
          legend.position = "none")+
    geom_point(data=lmdata, aes(x=PC1, y=MDS1, fill=population),
               col="black", size=7, shape=21) +
    scale_fill_manual(values=mycols) +
    xlab(xlab) + ylab(ylab) +
    geom_smooth(data=lmdata, aes(x=PC1, y=MDS1), method=lm) +
    annotate("text", x=(min(lmdata$PC1)+max(lmdata$PC1))/2,y=min(lmdata$MDS)*1.2, 
             label=lm_eqn(data.frame(x=lmdata$PC1, y=lmdata$MDS1)),
             parse = TRUE, size=5)+
    ggtitle(paste0(pop," (n=", nrow(lmdata),")"))
  
  print(out)
  
  len <- length(lm_plot_out)
  nm <- c(names(lm_plot_out),pop)
  lm_plot_out[[len+1]] <- out
  
  names(lm_plot_out) <- nm
  
  #fileout <- paste0("output_files/revised_MAG/FigureS3/FigS3D",pop)
  fileout="output_files/revised_MAG/revison_fig/FigureS3/FigS3D"
  
  for(d in c("svg","pdf","png")){
    ggsave(plot=out, filename = paste0(fileout,".",d),device = d,width = 5,height = 5,units = "in",scale = 1.5)
  }
  saveRDS(out, file = paste0(fileout,".rds"))
  
}


library(gridExtra)
col <- grid.arrange(lm_plot_out[[1]],
                    lm_plot_out[[2]],
                    lm_plot_out[[3]],
                    lm_plot_out[[4]],
                    layout_matrix=rbind(1:4))

#outdir <- paste0("output_files/revised_MAG/FigureS3",pop)

outdir="output_files/revised_MAG/revison_fig/FigureS3"
saveRDS(col, file.path(outdir,"FigS3D.rds"))
ggsave(plot=col, filename = file.path(outdir,"FigS3D.png"),height = 5,width = 19,device = "png", units = "in",scale = 0.9)
ggsave(plot=col, filename = file.path(outdir,"FigS3D.svg"),height = 5,width = 19,device = "svg", units = "in",scale = 0.9)
ggsave(plot=col, filename = file.path(outdir,"FigS3D.pdf"),height = 5,width = 19,device = "pdf", units = "in",scale = 0.9)


