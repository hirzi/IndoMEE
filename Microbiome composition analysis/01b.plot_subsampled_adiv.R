rm(list=ls())

# packages
library(ggplot2) 
library(ggpubr)
library(scales)
library(gridExtra)
library(reshape2)
library(dplyr)
source("99.community_data_utils.R") # load custom scripts

###########
## Paths ##
###########
setwd("C:/Users/caf77/OneDrive - University of Cambridge/Documents/Analysis/MOBILE/0_Pilot_reviewed")
figout_dir="C:/Users/caf77/OneDrive - University of Cambridge/Documents/Analysis/MOBILE/0_Pilot_reviewed/figout/adiv"
dir.create(figout_dir,recursive=T,showWarnings=F)

#################################
## Setups variables and levels ##
#################################
# load metadata
source("99.metadata_with_SPMP.R") 
head(mypop)
print(pop.ord)
print(mycols)
print(myshape)

# load dataset
df0 <- read.delim("input_files/alpha_diversity_subsampled_wFD.tsv")
colnames(df0)[6]<-"simpson_E"

# samples that did not pass sub-sampling
exclude.list <- scan("input_files/adiv_excluded_samples.txt",what=character())

# set limits
metric_max <- apply(df0[2:8],2,
                    function(x) {
                      y=max(x)*1.1
                      ifelse(y<1, round(y,1), 
                             ifelse(y<10,round(y,0),ceiling(y/10)*10))
                    })

metric_min <- c(0,0.4,0.4,0,0.02,0.04,0)
names(metric_min) <- names(metric_max)

metric_step <- unlist(lapply(metric_max, function(x){
ifelse(x<1,0.1,
       ifelse(x<2, 0.25,
         ifelse(x<5,0.5,
                ifelse(x<10,2,
                       ifelse(x<150,20,100)))))}))
metric_step["simpson"]  <- 0.1

# NA check
df0[is.na(df0$observed),]
df0[is.na(df0$shannon),]
df0[is.na(df0$simpson),]
df0[is.na(df0$simpsonE),]
df0[is.na(df0$dominance),]
df0[is.na(df0$faith),]

# Plot and Save
iter_list <- unique(df0$subsample)
for(iset in iter_list){
  print(iset)
  adiv.df <- subset(df0, df0$subsample==iset)
  adiv.df <- merge.data.frame(mypop[,c("sampleid","population","urbanisation","n","group")], adiv.df[,-length(adiv.df)], by="sampleid")                      
  adiv.df.long <- reshape2::melt(adiv.df, id.vars=c("sampleid","population","urbanisation","n","group"), 
                                 variable.name="metrics",value.name="value")
  out_dir <- paste0(figout_dir,"/",iset)
  dir.create(out_dir,recursive=T,showWarnings=F)

  for(i in unique(adiv.df.long$metrics)){
    print(i)
    sumdata <- subset(adiv.df.long, metrics==i)
    not_compared <- subset(sumdata, sampleid%in%exclude.list)
    sumdata <- subset(sumdata, !sampleid%in%exclude.list)
    
    # set factors and sort
    sumdata$population <- factor(sumdata$population, pop.ord)
    not_compared$population <- factor(not_compared$population, pop.ord)
    sumdata <- sumdata[order(sumdata$population),]
    not_compared <- not_compared[order(not_compared$population),]
    
    # set width dodge
    not_compared$dodge <- ifelse(not_compared$population=="Balinese_DPS",-0.23,
                                 ifelse(not_compared$population=="Indian_SIN", 0.08,
                                        ifelse(not_compared$population=="Chinese_SIN",0.22,0)))
    not_compared$xpos <- not_compared$dodge+as.numeric(not_compared$urbanisation)
    
    # compare stats
    urbpair <- list(c("rural", "remote"),
                    c("urban","rural"),
                    c("urban","remote"))
    
    mytest <- pairwise.wilcox.test(sumdata$value,sumdata$urbanisation,p.adjust.method="fdr")$p.value
    
    mypval <- c(mytest[urbpair[[1]][1],urbpair[[1]][2]],
                mytest[urbpair[[2]][1],urbpair[[2]][2]],
                mytest[urbpair[[3]][1],urbpair[[3]][2]])
    
    nSignif <- mypval<0.05
    my.annot <-  list(pairs=urbpair[nSignif],
                      p.value=mypval[nSignif],
                      p.stars=gtools::stars.pval(mypval)[nSignif],
                      star.height=c(metric_max[i], 1.1*metric_max[i], 1.2*metric_max[i])[nSignif]
                      )
    
    lim=metric_max[i]*(1.3)
    ybreak=seq(0,metric_max[i],by=metric_step[i])
    
    p <- ggplot(sumdata,aes(x=urbanisation, y=value)) + 
      theme_bw(base_size=20,base_rect_size=2, base_line_size=2) +
      theme(panel.grid=element_blank(), title=element_text(size=14),
            legend.key.size=unit(7,"mm"), legend.text=element_text(size=14), legend.title=element_text(size=16))+
      geom_boxplot(size=1, col="grey50", outliers=F) +
      geom_signif(comparisons=my.annot$pairs, map_signif_level=T,
                  y_position=my.annot$star.height,
                  annotations=my.annot$p.stars,
                  size=1,textsize=6)+
      geom_point(data=sumdata, aes(x=urbanisation, y=value, colour=population, fill=population, shape=urbanisation, size=urbanisation), 
                 position=position_jitterdodge(jitter.width=0.2, dodge.width=0.6, seed=1))+
      geom_point(data=not_compared, aes(x=xpos, y=value, group=population), 
                 fill="white", colour="black", shape=myshape[not_compared$urbanisation], size=3)+
      scale_y_continuous(expand=c(0.07,0), 
                         limits=c(metric_min[i], lim),
                         breaks=ybreak) +
      scale_color_manual(values=mycols) +
      scale_fill_manual(values=mycols) +
      scale_shape_manual(values=myshape) +
      scale_size_manual(values=c(3,3,2.7)) +
      ylab(i) +
      ggtitle(iset)
    
    # save
    filepath <- paste0(out_dir,"/", i, ".svg")
    ggsave(filename=filepath, plot=p, device="svg",width=8, height=6,units="in", dpi="print", bg="white")
    
    filepath <- paste0(out_dir,"/", i, ".png")
    ggsave(filename=filepath, plot=p, device="png",width=8, height=6,units="in", dpi=600, bg="white")
    
    filepath <- paste0(out_dir,"/", i, ".rds")
    saveRDS(p + theme(axis.title=element_blank(),
                      axis.text=element_text(size=15), 
                      line=element_line(size=1), 
                      panel.border=element_rect(linewidth=1.2)), filepath)
    #print(p)

  }
}

# collate plot
out_dir=paste0(figout_dir,"/collated")
dir.create(out_dir,recursive=T,showWarnings=F)

folder_list <- list.files(figout_dir,pattern="hybrid|comboRef")
metric_list <- as.character(unique(adiv.df.long$metrics))
for(i in metric_list){
  print(i)
  myplot <- lapply(folder_list, function(x) {
    p<-paste0(paste(figout_dir,x,i,sep="/"),".rds")
    readRDS(p)
    })
  myplot <- myplot[c(3,4,1,2)]
  p_col <- annotate_figure(p=ggarrange(plotlist=myplot, ncol=2, nrow=2,
                                         labels="AUTO", font.label=list(size=20),
                                         common.legend=T,legend="right"),
                           top=text_grob(toupper(i), face="bold", size=20),
                           bottom=text_grob("Subsampling was set at an even depth of 29.8 million reads, and samples with fewer reads were excluded (2 DPS, 6 SIN).", 
                                            color="grey55", hjust=0.5, x=0.5, size=11))
  
  filepath <- paste0(out_dir,"/", i, ".svg")
  ggsave(filename=filepath, plot=p_col, device="svg",width=6, height=6,units="in", dpi="print", bg="white", scale=1.8)
  
  filepath <- paste0(out_dir,"/", i, ".png")
  ggsave(filename=filepath, plot=p_col, device="png",width=6, height=6,units="in", dpi=600, bg="white", scale=1.8)
  
  #print(p_col)
}


# by population
legend_theme <- theme(
  legend.background=element_blank(),
  legend.margin=margin(0, 0, 0, 0),      # inner margin inside legend box
  legend.box.margin=margin(0, 0, 0, 0),  # space around the legend box
  legend.spacing=unit(5, "pt"),          # spacing between legend items
  legend.key.height=unit(4, "mm"),       # shrink legend key if needed
  legend.key.width=unit(4, "mm"),
  legend.key.size=unit(3,"mm"), 
  legend.text=element_text(size=11), 
  legend.title=element_text(size=14),
  legend.position="right",              # or "top"
  legend.justification="centre"         # aligns to top of the margin
)

log10rev <- trans_new(
  name="log10rev",
  transform=function(x) -log10(x),
  inverse=function(x) 10^(-x),
  breaks=function(x) {
    pbreaks <- c(1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12)
    pbreaks <- pbreaks[pbreaks >= min(10^(-x)) & pbreaks <= max(10^(-x))]
    -log10(pbreaks)
  },
  format=function(x) format(10^(-x), scientific=TRUE)
)

iter_list <- unique(df0$subsample)
for(iset in iter_list){
  print(iset)
  adiv.df <- subset(df0, df0$subsample==iset)
  adiv.df <- merge.data.frame(mypop[,c("sampleid","population","urbanisation","n","group")], adiv.df[,-length(adiv.df)], by="sampleid")                      
  adiv.df.long <- reshape2::melt(adiv.df, id.vars=c("sampleid","population","urbanisation","n","group"), 
                                 variable.name="metrics",value.name="value")
  out_dir <- paste0(figout_dir,"/",iset)
  dir.create(out_dir,recursive=T,showWarnings=F)
  
  for(i in unique(adiv.df.long$metrics)){
    print(i)
    sumdata <- subset(adiv.df.long, metrics==i)
    sumdata$population <- factor(sumdata$population, pop.ord)
    
    # exclude samples from main dataset
    sumdata <- subset(adiv.df.long, metrics==i)
    not_compared <- subset(sumdata, sampleid%in%exclude.list)
    sumdata <- subset(sumdata, !sampleid%in%exclude.list)
  
    # set factors and sort
    sumdata$population <- factor(sumdata$population, pop.ord)
    not_compared$population <- factor(not_compared$population, pop.ord)
    sumdata <- sumdata[order(sumdata$population),]
    not_compared <- not_compared[order(not_compared$population),]
    
    # set limits and breaks
    height_step=max(sumdata$value)*0.03
    lim=metric_max[i]+(height_step*length(combn(pop.ord,2))*0.7)
    ybreak=seq(0,metric_max[i],by=metric_step[i])
    max_width <- max(nchar(pop.ord))
    padded.ybreak<-sprintf(paste0("%", max_width+5, "s"),as.character(ybreak))
    
    p1 <-  ggplot(sumdata,aes(x=population, y=value)) + 
      theme_bw(base_size=14,base_rect_size=1, base_line_size=1) +
      theme(plot.margin=margin(t=0, r=0, b=0, l=0, unit="pt"),
            panel.background=element_blank(),
            panel.grid=element_blank(), 
            panel.grid.major.x=element_line(size=0.8, colour="grey80", linetype=3),
            axis.title=element_blank(),
            axis.text.y=element_text(size=12),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            title=element_text(size=14)) +
      legend_theme +
      geom_boxplot(size=1, col="grey50", outliers=F) +
      geom_point(data=sumdata, aes(x=population, y=value, colour=population, fill=population, shape=urbanisation, size=urbanisation), 
                 position=position_jitterdodge(jitter.width=0.2, dodge.width=0.6, seed=1))+
      geom_point(data=not_compared,aes(x=population, y=value, shape=urbanisation), fill="white",colour="black", size=3)+
      scale_y_continuous(expand=c(0.07,0),
                         limits=c(metric_min[i],metric_max[i]),
                         breaks=ybreak,labels=padded.ybreak) +
      scale_color_manual(values=mycols) +
      scale_fill_manual(values=mycols) +
      scale_shape_manual(values=myshape) +
      scale_size_manual(values=c(3,3,2.7)) +
      ylab(i) +
      ggtitle(paste0(iset," - ",i))
    
    # get significance
    mytest <- pairwise.wilcox.test(sumdata$value,sumdata$population,p.adjust.method="fdr",paired=F, exact=F)$p.value
    pval.df <- melt(mytest)
    pval.df <- na.omit(pval.df)
    mirror <- pval.df
    colnames(mirror)[c(1,2)] <- c("Var2","Var1")
    mirror <- mirror[,c(2,1,3)]
    
    pval.df <- unique(rbind(pval.df,
                            mirror))
    
    pval.df$log_pval <- -log(pval.df$value,base=10)
    pval.df$Var1 <- factor(pval.df$Var1, rev(pop.ord))
    pval.df$Var2 <- factor(pval.df$Var2, pop.ord)
    
    p2 <- ggplot(pval.df)+
      theme_bw(base_size=14,base_rect_size=1, base_line_size=1) +
      theme(plot.margin=margin(t=0, r=0, b=0, l=0, unit="pt"),
            panel.background=element_blank(),
            panel.grid=element_blank(),
            axis.title=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_text(size=10),
            axis.ticks.x=element_blank(),
            title=element_text(size=14),
            theme(panel.spacing=unit(0.5, "cm")))+
      legend_theme+
      geom_tile(aes(x=Var2,y=Var1,fill=value))+
      scale_fill_gradient2(
        low="white",
        mid="yellow",
        high="tomato",
        midpoint=0.05,
        trans=log10rev,
        breaks=c(1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12),
        labels=scales::trans_format("log10", scales::math_format(10^.x)),
        name="p-value"
      )
    p2
    
    # Extract both legends
    legend1 <- get_legend(p1)
    legend2 <- get_legend(p2)
    
    # Remove legends from plots
    p1_clean <- p1 + theme(legend.position="none")
    p2_clean <- p2 + theme(legend.position="none")
    
    # Combine legends vertically (can also do horizontally)
    both_legends <- ggarrange(legend1, legend2, ncol=1)

    
    # PNG and SVG
    p <- ggarrange(
      ggarrange(p1_clean, 
                p2_clean, 
                ncol=1, align="v",heights=c(3,2)),
      both_legends,
      ncol=2,
      widths=c(4, 1))
    
        # save
        filepath <- paste0(out_dir,"/", i, "_bypop.svg")
        ggsave(filename=filepath, plot=p, device="svg",width=8, height=6,units="in", dpi="print", bg="white")
        
        filepath <- paste0(out_dir,"/", i, "_bypop.png")
        ggsave(filename=filepath, plot=p, device="png",width=8, height=6,units="in", dpi=600, bg="white")
    
    # RDS
    p <- ggarrange(p1_clean+ggtitle(iset),
                   p2_clean,
                   ncol=1, 
                   align="v",
                   heights=c(3,2))
        # save
        filepath <- paste0(out_dir,"/", i, "_bypop.rds")
        saveRDS(p, filepath)
        dir.create(paste0(figout_dir,"/collated"),recursive=T,showWarnings=F)
        filepath <- paste0(figout_dir,"/collated/unified_legend.rds")
        saveRDS(both_legends, filepath)
    
    # print(p) # display output
    }
    }

# collate plot
out_dir=paste0(figout_dir,"/collated")
dir.create(out_dir,recursive=T,showWarnings=F)

folder_list <- list.files(figout_dir,pattern="hybrid|comboRef")
folder_list <- folder_list[c(3,4,1,2)]
metric_list <- as.character(unique(adiv.df.long$metrics))

for(i in metric_list){
  print(i)
  myplot <- lapply(folder_list, function(x) {
    p<-paste0(paste(figout_dir,x,i,sep="/"),"_bypop.rds")
    readRDS(p)
  })
  
  p_col <- annotate_figure(
    ggarrange(
    ggarrange(plotlist=myplot, ncol=2, nrow=2,
              labels="AUTO", font.label=list(size=20),
              common.legend=T,legend="right"),
    both_legends,
    ncol=2,
    widths=c(4, 1)),
    top=text_grob(toupper(i), face="bold", size=20),
    bottom=text_grob("Subsampling was set at an even depth of 29.8 million reads, and samples with fewer reads were excluded (2 DPS, 6 SIN).", color="grey55", hjust=0.5, x=0.5, size=11))
  
  filepath <- paste0(out_dir,"/", i, "_bypop.svg")
  ggsave(filename=filepath, plot=p_col, device="svg",width=6, height=6.5,units="in", dpi="print", bg="white", scale=1.8)
  
  filepath <- paste0(out_dir,"/", i, "_bypop.png")
  ggsave(filename=filepath, plot=p_col, device="png",width=6, height=6.5,units="in", dpi=600, bg="white", scale=1.8)
  
  #print(p_col)
}


## Figure2A and Figure S2
select.metrics <- c("shannon","observed","faith")
for(i in select.metrics){
  print(i)
  sumdata <- subset(adiv.df.long, metrics==i)
  
  urbpair <- list(c("rural", "remote"),
                  c("urban","rural"),
                  c("urban","remote"))
  
  mytest <- pairwise.wilcox.test(sumdata$value,sumdata$urbanisation,p.adjust.method="fdr")$p.value
  
  mypval <- c(mytest[urbpair[[1]][1],urbpair[[1]][2]],
              mytest[urbpair[[2]][1],urbpair[[2]][2]],
              mytest[urbpair[[3]][1],urbpair[[3]][2]])
  
  my.annot <-  list(pairs=urbpair,
                    p.value=mypval,
                    p.stars=gtools::stars.pval(mypval),
                    star.height=c(1.1*max(sumdata$value), 1.18*max(sumdata$value), 1.26*max(sumdata$value)))
  
  p <- ggplot(sumdata,aes(x=urbanisation, y=value)) + 
    theme_bw(base_size=12,base_rect_size=2, base_line_size=0.8) +
    theme(panel.grid=element_blank())+
    geom_boxplot(size=0.8, col="grey50", outliers=F) +
    geom_signif(comparisons=my.annot$pairs[my.annot$p.value<0.05], 
                y_position=my.annot$star.height[my.annot$p.value<0.05],
                annotations=my.annot$p.stars[my.annot$p.value<0.05],
                size=0.8,textsize=4)+
    geom_point(data=sumdata, aes(x=urbanisation, y=value, colour=population, fill=population, shape=urbanisation, size=urbanisation), 
               position=position_jitterdodge(jitter.width=0.2, dodge.width=0.6, seed=1))+
    scale_y_continuous(expand=c(0.07,0)) +
    scale_color_manual(values=mycols) +
    scale_fill_manual(values=mycols) +
    scale_shape_manual(values=myshape) +
    scale_size_manual(values=c(2,2,1.5)) + 
    theme(legend.position="none")+
    theme(axis.title.x=element_blank())+
    ylab(i)
  
  
  # save
  if(i=="shannon"){
    filepath <- paste0("figout/Figure2/Fig2A_", i, ".svg")
    ggsave(filename=filepath, plot=p, device="svg",width=3, height=6,units="in", dpi="print", bg="white")
    
    filepath <- paste0("figout/Figure2/Fig2A_", i, ".png")
    ggsave(filename=filepath, plot=p, device="png",width=3, height=6,units="in", dpi=600, bg="white")
    
    filepath <- paste0("figout/Figure2/Fig2A_", i, ".rds")
    saveRDS(p, filepath)
  }else{
    filepath <- paste0("figout/FigureS2/FigS2A_", i, ".svg")
    ggsave(filename=filepath, plot=p, device="svg",width=3, height=6,units="in", dpi="print", bg="white")
    
    filepath <- paste0("figout/FigureS2/FigS2A_", i, ".png")
    ggsave(filename=filepath, plot=p, device="png",width=3, height=6,units="in", dpi=600, bg="white")
    
    filepath <- paste0("figout/FigureS2/FigS2A_", i, ".rds")
    saveRDS(p, filepath)
  }
  
  print(p)
}


####################
# Collocation test #
####################
df0 <- read.delim("input_files/alpha_diversity_subsampled_wFD.tsv")
colnames(df0)[6]<-"simpson_E"
adiv <- data.frame(subset(df0, subsample=="hybrid_subsampled_iter0"),row.names = NULL)
compare_df <- merge.data.frame(mypop[,c("sampleid","population","urbanisation","n","group")], adiv, by="sampleid")

# Wilcoxon / Kruskal
pairwise.wilcox.test(compare_df$shannon, g = compare_df$urbanisation, p.adjust.method = "none")
pairwise.wilcox.test(compare_df$shannon, g = compare_df$population, p.adjust.method = "none") -> pop_wilcox
round(as.matrix(pop_wilcox$p.value)[7:10,8:10],4)

tmp <- compare_df[compare_df$population%in%c("Balinese_DPS","Malay_SIN","Chinese_SIN","Indian_SIN"),c("shannon","population")]
wilcox.test(shannon ~ population, data = subset(tmp, tmp$population%in%c("Balinese_DPS","Chinese_SIN")), alternative="greater")
wilcox.test(shannon ~ population, data = subset(tmp, tmp$population%in%c("Balinese_DPS","Malay_SIN")), alternative="greater")
wilcox.test(shannon ~ population, data = subset(tmp, tmp$population%in%c("Balinese_DPS","Indian_SIN")), alternative="greater")
wilcox.test(tmp$shannon ~ ifelse(tmp$population == "Balinese_DPS","no","yes"), alternative = "greater")


# HUL vs RPN (collocation - Borneo)
wilcox.test(shannon ~ population, alternative="greater",
            data=subset(compare_df, compare_df$population %in% c("Punan_HUL","Punan_RPN")))
t.test(shannon ~ population,
            data=subset(compare_df, compare_df$population %in% c("Punan_HUL","Punan_RPN")))
ks.test(compare_df$shannon[compare_df$population=="Punan_HUL"], 
        compare_df$shannon[compare_df$population=="Punan_RPN"])

# SUL vs RPN (collocation - Borneo)
wilcox.test(shannon ~ population,
            data=subset(compare_df, compare_df$population %in% c("Punan_SUL","Punan_RPN")))
t.test(shannon ~ population,
       data=subset(compare_df, compare_df$population %in% c("Punan_SUL","Punan_RPN")))
ks.test(compare_df$shannon[compare_df$population=="Punan_SUL"], 
        compare_df$shannon[compare_df$population=="Punan_RPN"])

# PDW vs DPS (collocation - Bali)
wilcox.test(shannon ~ population, alternative ="greater",
            data=subset(compare_df, compare_df$population %in% c("Balinese_PDW","Balinese_DPS")))
t.test(shannon ~ population,alternative ="greater",
       data=subset(compare_df, compare_df$population %in% c("Balinese_PDW","Balinese_DPS")))
ks.test(compare_df$shannon[compare_df$population=="Balinese_PDW"], 
        compare_df$shannon[compare_df$population=="Balinese_DPS"])

# PDW vs RPN (both rural)
wilcox.test(shannon ~ population,
            data=subset(compare_df, compare_df$population %in% c("Balinese_PDW","Punan_RPN")))
ks.test(compare_df$shannon[compare_df$population=="Balinese_PDW"], 
        compare_df$shannon[compare_df$population=="Punan_RPN"])

#  PDW vs BSP 
wilcox.test(shannon ~ population,
            data=subset(compare_df, compare_df$population %in% c("Balinese_PDW","Basap_BRU")))
ks.test(compare_df$shannon[compare_df$population=="Balinese_PDW"], 
        compare_df$shannon[compare_df$population=="Basap_BRU"])

#  PDW vs RPN 
wilcox.test(shannon ~ population, 
            data=subset(compare_df, compare_df$population %in% c("Balinese_PDW","Punan_RPN")))
ks.test(compare_df$shannon[compare_df$population=="Balinese_PDW"], 
        compare_df$shannon[compare_df$population=="Punan_RPN"])

#  3pop SIN - Kurskal -?
sin_df <- subset(compare_df, compare_df$population %in% c("Chinese_SIN","Malay_SIN","Indian_SIN"))
pairwise.wilcox.test(sin_df$shannon, g=sin_df$population, p.adjust.method = "fdr")
