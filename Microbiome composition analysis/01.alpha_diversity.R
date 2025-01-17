# custom functions
convert_to_CPM <- function(comm.data){
  comm.df <- comm.data # abundance matrix with columns as sampleid
  comm.df <- comm.df*10^6}
prevalence_filter <- function(x, n) {
  x <- reads # matrix of abundance with sampleid as columns
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

############################
## Load setups and levels ##
############################
mypop <- readRDS("input_files/mypop_wSPMP.rds")
pop.ord <- levels(mypop$population)
mycols <- read.delim("input_files/MAGs/MAGs_pop_mycols.tsv", header=T)
mycols <- mycols[mycols$group %in% unique(mypop$population),]
col.ord <- mycols$group
mycols <- mycols$mycols
names(mycols) <- col.ord
mycols <- mycols[pop.ord]

myshape <- c(21,22,24)
names(myshape) <- c("remote","rural","urban")

#################
## Input Files ##
#################
# abundance data
df <- read.delim("input_files/MAGs/revised/Hybrid_wSPMP/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv")
sampleid <- readRDS("output_files/sampleid_n116.rds")
sporeids <- readRDS("output_files/sampleid_spore.rds")
reads <- df[,colnames(df) %in% c(sampleid,sporeids)]
rownames(reads) <- df$original_bin
comm.df <- convert_to_CPM(reads)
table(colSums(comm.df)==0)
dim(comm.df)


# load tree
library("ape")
library("phangorn")
library("picante")

tree <- read.tree("input_files/MAGs/gtdbtk.bac120.user_msa.fasta.treefile")
tree <- as.phylo(tree)
tree <- midpoint(tree)

# check if tree matches to otu table
table(rownames(comm.df) %in% tree$tip.label)
not.present.table <- rownames(comm.df)[!rownames(comm.df) %in% tree$tip.label]
table(tree$tip.label %in% rownames(comm.df))
not.present.tree <- tree$tip.label[!tree$tip.label %in% rownames(comm.df)]

# prune tree (only if the tree contain more OTUs than the OTU table)
# tree.pruned <- prune.missing(t(comm.df), phylo = tree)[[2]]

# remove otus not present in tree
comm.df.archaea <- comm.df[!rownames(comm.df) %in% tree$tip.label,]
comm.df.bacteria <- comm.df[rownames(comm.df) %in% tree$tip.label,]

####################
## Alpha Diversiy ##
####################
library(abdiv)
shannon <- as.numeric(apply(comm.df.bacteria,2, function(x) abdiv::shannon(x,base = 10)))
simpson <- as.numeric(apply(comm.df.bacteria,2, function(x) abdiv::simpson(x)))
invsimpson <- as.numeric(apply(comm.df.bacteria,2, function(x) abdiv::invsimpson(x)))
simpsonE <- as.numeric(apply(comm.df.bacteria,2, function(x) abdiv::simpson_e(x)))
dominance <- as.numeric(apply(comm.df.bacteria,2, function(x) abdiv::strong(x)))
observed <- vegan::specnumber(comm.df.bacteria, MARGIN = 2)
faith <- as.numeric(apply(comm.df.bacteria,2, function(x) abdiv::faith_pd(abdiv::match_to_tree(x, tree=tree), tree=tree)))

adiv.df <- data.frame(sampleid = colnames(comm.df.bacteria),
                      observed, shannon, simpson, invsimpson, simpsonE, faith, dominance)

scount <- data.frame(table(mypop$population), stringsAsFactors = F)

colnames(scount) <- c("population","n")
scount$group <- paste0(scount$population,"\n(n=", scount$n,")")
mypop_n <- merge.data.frame(mypop, scount, by="population")

adiv.df <- merge.data.frame(mypop_n, adiv.df, by="sampleid")                      

adiv.df$urbanisation <- ifelse(adiv.df$population %in% c("Asmat_DAI","Punan_PBS","Punan_SUL","Punan_HUL"), "remote",
                         ifelse(adiv.df$population %in% c("Punan_RPN","Basap_BRU","Balinese_PDW"), "rural",
                                ifelse(adiv.df$population %in% c("Balinese_DPS","Chinese_SIN","Malay_SIN","Indian_SIN"), "urban","")))


#########
#compare#
#########
adiv.metrics <- c("observed","faith","shannon","simpson","invsimpson","simpsonE","dominance")
wilcox.out.pop <- apply(adiv.df[,adiv.metrics], 2, function(x) pairwise.wilcox.test(x, g = adiv.df$population,p.adjust.method = "fdr", exact=F))

signif.wilcox.pop <- list()
for(metric in names(wilcox.out.pop)){
  udf <- wilcox.out.pop[[metric]]
  pairout <- na.omit(reshape2::melt(data.frame(var1=rownames(udf$p.value), udf$p.value), id.vars="var1", variable.name="var2",value.name = "p.adj"))
  pairout.signif <- pairout[pairout$p.adj<0.05,]
  pairout.signif[,c(1,2)] <- apply(pairout.signif[,c(1,2)], 2, function(x) as.character(x))
  pairout.signif$pair<-paste(pairout.signif$var1,pairout.signif$var2,sep="|")
  rownames(pairout.signif) <- NULL
  
  pairout.signif.df <- data.frame()
  for(i in 1:nrow(pairout.signif)){
      x <- pairout.signif[i,c(1:2)]
      tmp <- subset(adiv.df, adiv.df$population %in% as.character(c(x[[1]],x[[2]])))
      out <- wilcox.test(tmp[,metric] ~ tmp$population, exact=F)
      dfout <- data.frame(pop1=paste(x[[1]]),
                     pop2=paste(x[[2]]),
                     pair=paste(x[[1]],x[[2]],sep="|"),
                     statistic=as.numeric(out$statistic),
                     p.value=as.numeric(out$p.value))
      dfout <- merge.data.frame(dfout, pairout.signif[,c(3,4)], by="pair")
      pairout.signif.df <- rbind(pairout.signif.df, dfout)
    }
    
  pairout.signif.df[,1:3] <- data.frame(apply(pairout.signif.df[,c(1:3)], 2, function(x) as.character(x)),stringsAsFactors = F)
  # pairout.signif.df$statistic <- gsub(pattern = "c\\(V\\ =\\ ",replacement = "",x = gsub(pattern = "\\)", replacement = "",x = pairout.signif.df$statistic))
  # pairout.signif.df$statistic <- as.numeric(pairout.signif.df$statistic)
  # pairout.signif.df$p.adjust <- as.numeric(pairout.signif.df$p.adjust)
  
  # View(pairout.signif.df)
  saveRDS(pairout.signif.df, file = paste("output_files/adiv/",metric,"_PopPair_Utest.rds",sep = ""))
  write.table(pairout.signif.df, file = paste("output_files/adiv/",metric,"_PopPair_Utest.tsv",sep = ""),quote = F, row.names = F, sep = "\t")
  
  len <- length(signif.wilcox.pop)+1
  signif.wilcox.pop[len] <- list(pairout.signif.df)
  names(signif.wilcox.pop)[len] <- metric
}
signif.wilcox.pop


wilcox.out.ls <- apply(adiv.df[,adiv.metrics], 2, function(x) pairwise.wilcox.test(x, g = adiv.df$urbanisation,p.adjust.method = "fdr"))

na.omit(reshape2::melt(data.frame(var1=rownames(wilcox.out.ls$observed$p.value), wilcox.out.ls$observed$p.value), id.vars="var1", variable.name="var2")) -> pairout
pairout.signif <- pairout[pairout$value<0.05,]
pairout.signif <- apply(pairout.signif, 2, function(x) as.character(x))
rownames(pairout.signif) <- NULL

pairout.signif.df <- data.frame()
for(i in 1:nrow(pairout.signif)){
  x <- pairout.signif[i,c(1:2)]
  tmp <- subset(adiv.df, adiv.df$urbanisation %in% as.character(c(x[[1]],x[[2]])))
  out <- wilcox.test(tmp$observed, g=tmp$urbanisation)
  dfout <- cbind(pop1=paste(x[[1]]),
                 pop2=paste(x[[2]]),
                 data.frame(pair=paste(x[[1]],x[[2]],sep="|")),
                 data.frame(rbind(out), row.names = i)[,c("statistic","p.value")])
  pairout.signif.df <- rbind(pairout.signif.df, dfout)
}
colnames(pairout.signif.df)[5] <- "p.adjust"
View(pairout.signif.df)

############
### PLOT ###
############
tmp <- adiv.df[,c("sampleid","population", "group", "urbanisation",adiv.metrics)]
adiv.df.long <- reshape2::melt(tmp, id.vars=c("sampleid","population", "group", "urbanisation"), variable.name="metrics",value.name="value")                      

library(ggplot2) 
library(dplyr)
for(i in unique(adiv.df.long$metrics)){
  print(i)
  sumdata <- subset(adiv.df.long, metrics==i)
  
  sumdata %>% group_by(population, group, urbanisation) %>% summarise(mean=mean(value),
                                            median = median(value),
                                            min = min(value),
                                            max = max(value),
                                            se=sd(value)/sqrt(length(value))) -> ggdata
  ggdata$lower <- ggdata$mean - ggdata$se
  ggdata$upper <- ggdata$mean + ggdata$se
  
  ggdata$urbanisation <- factor(ggdata$urbanisation, c("remote","rural","urban"))
  p <- ggplot(ggdata) + theme_bw(base_size = 20) +
    geom_point(aes(x=urbanisation, y=mean, fill=population, color=population, shape=urbanisation), size=3, stat="identity", position = position_dodge(width = 0.5))+
    geom_errorbar(aes(x=urbanisation, ymin=lower, ymax=upper, colour=population), width=0.2, position = position_dodge(width = 0.5)) +
    ggtitle(i) +
    scale_y_continuous(expand = c(0,0),limits = c(0,1.1*max(ggdata$upper))) +
    scale_shape_manual(values=myshape)+
    scale_color_manual(values= mycols)+
    scale_fill_manual(values= mycols)
  
  print(p)
  
  filepath <- paste0("output_files/revised_MAG/adiv/", i, ".svg")
  ggsave(filename = filepath, plot = p, device = "svg",width = 6, height = 6,units = "in", dpi="print", bg = "white")
  
  filepath <- paste0("output_files/revised_MAG/adiv/", i, ".png")
  ggsave(filename = filepath, plot = p, device = "png",width = 6, height = 6,units = "in", dpi=72, bg = "white")
  
}

library(ggpubr)
for(i in unique(adiv.df.long$metrics)){
  print(i)
  sumdata <- subset(adiv.df.long, metrics==i)
  
  urbpair <- list(c("rural", "remote"),
                  c("urban","rural"),
                  c("urban","remote"))

  mytest <- pairwise.wilcox.test(sumdata$value,sumdata$urbanisation,p.adjust.method = "fdr")$p.value
    
  mypval <- c(mytest[urbpair[[1]][1],urbpair[[1]][2]],
              mytest[urbpair[[2]][1],urbpair[[2]][2]],
              mytest[urbpair[[3]][1],urbpair[[3]][2]])
  
  my.annot <-  list(pairs=urbpair,
                   p.value=mypval,
                   p.stars=gtools::stars.pval(mypval),
                   star.height=c(1.1*max(sumdata$value), 1.18*max(sumdata$value), 1.26*max(sumdata$value)))

  p <- ggplot(sumdata,aes(x=urbanisation, y=value)) + 
    theme_bw(base_size = 20,base_rect_size = 2, base_line_size = 2) +
    theme(panel.grid=element_blank())+
    geom_boxplot(size=1, col="grey50", outliers = F) +
    geom_signif(comparisons = my.annot$pairs[my.annot$p.value<0.05], 
                y_position = my.annot$star.height[my.annot$p.value<0.05],
                annotations = my.annot$p.stars[my.annot$p.value<0.05],
                size=1.2,textsize = 6)+
    geom_point(data=sumdata, aes(x=urbanisation, y=value, colour=population, fill=population, shape=urbanisation, size=urbanisation), 
               position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6, seed=1))+
    scale_y_continuous(expand = c(0.07,0)) +
    scale_color_manual(values = mycols) +
    scale_fill_manual(values = mycols) +
    scale_shape_manual(values = myshape) +
    scale_size_manual(values = c(3,3,2.7)) +
    ylab(i)
  
  
  # save
  filepath <- paste0("output_files/revised_MAG/adiv/", i, ".svg")
  ggsave(filename = filepath, plot = p, device = "svg",width = 8, height = 6,units = "in", dpi="print", bg = "white")
  
  filepath <- paste0("output_files/revised_MAG/adiv/", i, ".png")
  ggsave(filename = filepath, plot = p, device = "png",width = 8, height = 6,units = "in", dpi=600, bg = "white")
  
  print(p)
  }

## Figure2A and Figure S2
select.metrics <- c("shannon","observed","faith")
for(i in select.metrics){
  print(i)
  sumdata <- subset(adiv.df.long, metrics==i)
  
  urbpair <- list(c("rural", "remote"),
                  c("urban","rural"),
                  c("urban","remote"))
  
  mytest <- pairwise.wilcox.test(sumdata$value,sumdata$urbanisation,p.adjust.method = "fdr")$p.value
  
  mypval <- c(mytest[urbpair[[1]][1],urbpair[[1]][2]],
              mytest[urbpair[[2]][1],urbpair[[2]][2]],
              mytest[urbpair[[3]][1],urbpair[[3]][2]])
  
  my.annot <-  list(pairs=urbpair,
                    p.value=mypval,
                    p.stars=gtools::stars.pval(mypval),
                    star.height=c(1.1*max(sumdata$value), 1.18*max(sumdata$value), 1.26*max(sumdata$value)))
  
  p <- ggplot(sumdata,aes(x=urbanisation, y=value)) + 
    theme_bw(base_size = 12,base_rect_size = 2, base_line_size = 0.8) +
    theme(panel.grid=element_blank())+
    geom_boxplot(size=0.8, col="grey50", outliers = F) +
    geom_signif(comparisons = my.annot$pairs[my.annot$p.value<0.05], 
                y_position = my.annot$star.height[my.annot$p.value<0.05],
                annotations = my.annot$p.stars[my.annot$p.value<0.05],
                size=0.8,textsize = 4)+
    geom_point(data=sumdata, aes(x=urbanisation, y=value, colour=population, fill=population, shape=urbanisation, size=urbanisation), 
               position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6, seed=1))+
    scale_y_continuous(expand = c(0.07,0)) +
    scale_color_manual(values = mycols) +
    scale_fill_manual(values = mycols) +
    scale_shape_manual(values = myshape) +
    scale_size_manual(values = c(2,2,1.5)) + 
    theme(legend.position = "none")+
    theme(axis.title.x = element_blank())+
    ylab(i)
  
  
  # save
  if(i=="shannon"){
    filepath <- paste0("output_files/revised_MAG/Figure2/Fig2A_", i, ".svg")
    ggsave(filename = filepath, plot = p, device = "svg",width = 3, height = 6,units = "in", dpi="print", bg = "white")
    
    filepath <- paste0("output_files/revised_MAG/Figure2/Fig2A_", i, ".png")
    ggsave(filename = filepath, plot = p, device = "png",width = 3, height = 6,units = "in", dpi=600, bg = "white")
    
    filepath <- paste0("output_files/revised_MAG/Figure2/Fig2A_", i, ".rds")
    saveRDS(p, filepath)
  }else{
  filepath <- paste0("output_files/revised_MAG/FigureS2/FigS2A_", i, ".svg")
  ggsave(filename = filepath, plot = p, device = "svg",width = 3, height = 6,units = "in", dpi="print", bg = "white")
  
  filepath <- paste0("output_files/revised_MAG/FigureS2/FigS2A_", i, ".png")
  ggsave(filename = filepath, plot = p, device = "png",width = 3, height = 6,units = "in", dpi=600, bg = "white")
  
  filepath <- paste0("output_files/revised_MAG/FigureS2/FigS2A_", i, ".rds")
  saveRDS(p, filepath)
  }
  
  print(p)
}
