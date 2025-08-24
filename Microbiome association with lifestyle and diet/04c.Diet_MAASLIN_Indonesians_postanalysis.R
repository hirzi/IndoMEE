library(ggplot2)
library(ggpubr)
library(ggvenn)

# load data
input.comm <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_comm.tsv")
input.metadata <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_metadata.tsv")
maaslin.taxonomy <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_full_taxonomy.tsv")
maaslin.taxonomy$Species_label <- gsub("_species","_sp.", maaslin.taxonomy$Species_deduplicated) # shorthand taxonomy

####################
## Venn Analysis ###
####################
library(ggvenn)
diet.cols <- c("#2E8B57","#DB7093","#63B8FF")
#names(diet.cols) <- c("Sago","Pork","Chicken")

### All data (n=73)
n <- length(unique(input.metadata$sampleid))
multivar.df <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/Multivar_diet+age+sex+bmi/all_results.tsv")
multivar.df <- multivar.df[multivar.df$qval<0.05,]

vs <- list(Sago=multivar.df$feature[multivar.df$value=="Sago"],
           Pork=multivar.df$feature[multivar.df$value=="Pork"],
           Chicken=multivar.df$feature[multivar.df$value=="Chicken"])

p <- ggvenn(vs,text_size = 4, set_name_size = 5, fill_color = diet.cols) + 
  labs(title = "Multivariate - All", subtitle=paste("Sample Size:",n), caption = "Adjusment: Age, Sex, BMI")
print(p)

ggsave(p, file="figout/FigureS3/FigS3F_1.svg",device = "svg",width = 8,height = 8.2, scale = 0.4)
ggsave(p, file="figout/FigureS3/FigS3F_1.pdf",device = "pdf",width = 8,height = 8.2, scale = 0.4)
ggsave(p, file="figout/FigureS3/FigS3F_1.png",device = "png",width = 8,height = 8.2, scale = 0.4)
saveRDS(p, file="figout/FigureS3/FigS3F_1.rds")

### Non-Bali
n <- length(unique(input.metadata$sampleid[!input.metadata$population %in% c("Balinese_DPS","Balinese_PDW")]))
multivar.df <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/Multivar_diet+age+sex+bmi_nonBali/all_results.tsv")
multivar.df <- multivar.df[multivar.df$qval<0.05,]

vs <- list(Sago=multivar.df$feature[multivar.df$value=="Sago"],
           Pork=multivar.df$feature[multivar.df$value=="Pork"],
           Chicken=multivar.df$feature[multivar.df$value=="Chicken"])

p <- ggvenn(vs,text_size = 4, set_name_size = 5, fill_color = diet.cols) + 
  labs(title = "Multivariate - DAI,BRU,HUL", subtitle=paste("Sample Size:",n), caption = "Adjusment: Age, Sex, BMI")
print(p)


ggsave(p, file="figout/FigureS3/FigS3F_2.svg",device = "svg",width = 8,height = 8.2, scale = 0.4)
ggsave(p, file="figout/FigureS3/FigS3F_2.pdf",device = "pdf",width = 8,height = 8.2, scale = 0.4)
ggsave(p, file="figout/FigureS3/FigS3F_2.png",device = "png",width = 8,height = 8.2, scale = 0.4)
saveRDS(p, file="figout/FigureS3/FigS3F_2.rds")

### population vs diet
#### All Samples (n=73)
n <- length(unique(input.metadata$sampleid))
base_dir <- "C:/Users/caf77/OneDrive - University of Cambridge/Documents/Analysis/MOBILE/0_Pilot_reviewed/output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/"
tests_list <- list.files(base_dir)
se <- c(grep("Multi",tests_list),
       grep("population",tests_list))
selected_tests <- unique(tests_list[se])
print(selected_tests)

selected_pair <- c("population_noranef","Multivar_diet+age+sex+bmi")
vs2 <- lapply(selected_pair, function(x) {
  x <- read.delim(paste0(base_dir,x,"/all_results.tsv"))
  x <- x[x$qval < 0.05,]
  x <- unique(x$feature[x$metadata %in% c("Sago","Pork","Chicken","population")])
  return(x)
  })

names(vs2) <- c("Population","Diet")

p2 <- ggvenn(vs2,text_size = 4,set_name_size = 5, fill_color = c("#7137c880","#FA8072")) +
  labs(title = "Multivariate Diet - All",subtitle=paste("Sample Size:",n), caption = "Adjusment: Age, Sex, BMI")
print(p2)

ggsave(p2, file="figout/Figure3/Fig3D.svg",device = "svg",width = 8,height = 8, scale = 0.4)
ggsave(p2, file="figout/Figure3/Fig3D.pdf",device = "pdf",width = 8,height = 8, scale = 0.4)
ggsave(p2, file="figout/Figure3/Fig3D.png",device = "png",width = 8,height = 8, scale = 0.4)
saveRDS(p2, file="figout/Figure3/Fig3D.rds")

#### Non-Bali only
n <- length(unique(input.metadata$sampleid[!input.metadata$population %in% c("Balinese_DPS","Balinese_PDW")]))
selected_pair <- c("population_noranef_nonBali_adj","Multivar_diet+age+sex+bmi_nonBali")
vs2 <- lapply(selected_pair, function(x) {
  x <- read.delim(paste0(base_dir,x,"/all_results.tsv"))
  x <- x[x$qval < 0.05,]
  x <- unique(x$feature[x$metadata %in% c("Sago","Pork","Chicken","population")])
  return(x)
})

names(vs2) <- c("Population","Diet")

p2 <- ggvenn(vs2,text_size = 4,set_name_size = 5, fill_color = c("#7137c880","#FA8072")) +
  labs(title = "Multivariate Diet - ASM,BRU,HUL", subtitle=paste("Sample Size:",n), caption = "Adjusment: Age, Sex, BMI")
print(p2)

ggsave(p2, file="figout/FigureS3/FigS3F_4.svg",device = "svg",width = 8,height = 8, scale = 0.4)
ggsave(p2, file="figout/FigureS3/FigS3F_4.pdf",device = "pdf",width = 8,height = 8, scale = 0.4)
ggsave(p2, file="figout/FigureS3/FigS3F_4.png",device = "png",width = 8,height = 8, scale = 0.4)
saveRDS(p2, file="figout/FigureS3/FigS3F_4.rds")

#############
## Barplot ##
#############
library(purrr)
base_dir <- "C:/Users/caf77/OneDrive - University of Cambridge/Documents/Analysis/MOBILE/0_Pilot_reviewed/output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/"

count_data <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/0_SUMMARY/count_signif_multivar_wnonBali.tsv")
count_data <- count_data[!count_data$model %in% c("old_base","old_population","Multivar_diet+age+sex"),]
count_data <- count_data[-grep("year|wConf", count_data$model),]

# A) Map variables
mapvar <- data.frame(model=count_data$model,
           Diet=ifelse(count_data$var=="population",0,1),
           #Sago=ifelse(count_data$var=="population",0,1),
           #Pork=ifelse(count_data$var=="population",0,1),
           #Chicken=ifelse(count_data$var=="population",0,1),
           Population=ifelse(count_data$var=="population",1,0),
           #Asmat_DAI=ifelse(count_data$var=="population",1,0),
           #Punan_HUL=ifelse(count_data$var=="population",1,0),
           #Balinese_DPS=ifelse(count_data$var=="population",1,0),
           #Balinese_PDW=ifelse(count_data$var=="population",1,0),
           Age_Sex_BMI=ifelse(count_data$adjustment!="none",1,0),
           #Age=ifelse(count_data$adjustment!="none",1,0),
           #Sex=ifelse(count_data$adjustment!="none",1,0),
           #BMI=ifelse(count_data$adjustment!="none",1,0),
           Population_ranef=ifelse(count_data$ranef=="none",0,1),
           Include_Balinese=as.numeric(unlist(lapply(count_data$model, function(x) {ifelse(length(grep("nonBali",x))>0, 0, 1)}))))
#mapvar$Balinese_DPS[mapvar$Include_Balinese==0] <- 0
#mapvar$Balinese_PDW[mapvar$Include_Balinese==0] <- 0
mapvar[mapvar$model=="Multivar_diet+age+sex+bmi+population","Diet"] <- 1
mapvar[mapvar$model=="Multivar_diet+age+sex+bmi+population_nonBali","Diet"] <- 1
#mapvar[mapvar$model=="Multivar_diet+age+sex+bmi+population",c("Sago","Pork","Chicken")] <- 1
#mapvar[mapvar$model=="Multivar_diet+age+sex+bmi+population_nonBali",c("Sago","Pork","Chicken")] <- 1

#mapvar <- mapvar[order(mapvar$Age),]
#mapvar <- mapvar[order(mapvar$Asmat_DAI),]
#mapvar <- mapvar[order(mapvar$Sago, decreasing = T),]
mapvar <- mapvar[order(mapvar$Age_Sex_BMI),]
mapvar <- mapvar[order(mapvar$Diet, decreasing = T),]
mapvar <- mapvar[order(mapvar$Population),]
mapvar <- mapvar[order(mapvar$Population_ranef),]
mapvar <- mapvar[order(mapvar$Include_Balinese, decreasing = T),]
rownames(mapvar) <- NULL
mapvar <- mapvar[,c("model","Diet","Population","Population_ranef","Age_Sex_BMI", "Include_Balinese")]



mapvar.long <- reshape2::melt(mapvar, id.vars="model")
mapvar.long$variable <- factor(mapvar.long$variable, rev(unique(mapvar.long$variable)))
mapvar.long$model <- factor(mapvar.long$model, mapvar$model)

# annotate signif
tmp <- mapvar.long
tmp$signif <- 0
for(mod in unique(count_data$model)){
  infile <- paste0(base_dir,mod,"/all_results.tsv")
  infile <- read.delim(infile)
  infile <- infile[infile$qval<0.05,]
  var <- unique(infile$value)
  rm(infile)
  if(!is.null(var)){
    var[grep(pattern = "Sago|Pork|Chicken", var)] <- "Diet"
    var[grep(pattern = "HUL|ASM|DPS|PDW", var)] <- "Population"
    var[grep(pattern = "bmi", var)] <- "Age_Sex_BMI"
    var[grep(pattern = "M", var)] <- "Age_Sex_BMI"
    var[grep(pattern = "Age", var)] <- "Age_Sex_BMI"
    var <- unique(var)
  } else {
    print("Error: no significant variables found. Please heck input file.")
    break
  }
  
  isPopRanef <- grep("popranef",mod)
  
  if(length(isPopRanef) == 0){
    for(i in var){tmp[tmp$model==mod&tmp$variable==i,"signif"] <- 1}
  } else {
    for(i in var){tmp[tmp$model==mod&tmp$variable==i,"signif"] <- 1}
    tmp[tmp$model==mod&tmp$variable=="Population_ranef","signif"] <- 1
  }
 
}

a <- ggplot(mapvar.long, aes(x=model)) + theme_minimal()+
  geom_dotplot(aes(y=variable, colour=factor(value), fill=factor(value)),dotsize = 0.6)+
  scale_fill_manual(values = c("grey90","black"))+
  scale_colour_manual(values = c("grey90","black"))+
  scale_y_discrete(expand = c(0,0), labels=rev(c("Diet","Population (fixed)", "Population (ranef)", "Age+Sex+BMI","Dataset:Include_Balinese")))+
  theme(panel.grid = element_blank())+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title=element_blank()) +
  theme(axis.text.y = element_text(vjust = -0.4, size=10))+
  theme(legend.position = "none")
a

ggplot(tmp, aes(x=model)) + theme_minimal()+
  geom_point(aes(y=variable, colour=factor(signif), fill=factor(value), shape=factor(signif)), size=3)+
  scale_shape_manual(values = c(16,21))+
  scale_fill_manual(values = c("grey90","black"))+
  scale_colour_manual(values = c("grey90","black"))+
  scale_y_discrete(expand = c(0,0), labels=rev(c("Diet","Population (fixed)", "Population (ranef)", "Age+Sex+BMI","Dataset:Include_Balinese")))+
  theme(panel.grid = element_blank())+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title=element_blank()) +
  theme(axis.text.y = element_text(vjust = -0.4, size=10))+
  theme(legend.position = "none")

# B) Barplot of feature count
multi_list <- lapply(paste0(base_dir,mapvar$model,"/all_results.tsv"), function(x) read.delim(x))
names(multi_list) <- mapvar$model
multi_list <- lapply(multi_list, function(x) {
  x <- x[x$qval < 0.05,]
  x$metadata[x$value %in% c("Sago","Pork","Chicken")] <- "Sago+Pork+Chicken"
  x$metadata[x$valur %in% c("Age","Gender","BMI")] <- "Age+Sex+BMI"
  return(x)
})

summarize_feature_by_value <- function(df) {
  value_counts <- table(df$value)
  features_by_value <- split(df$feature, df$value)
  
  # Initialize empty outputs for shared feature stats
  shared_counts <- numeric(0)
  shared_lists <- list()
  
  if (length(features_by_value) >= 2) {
    combos <- combn(names(features_by_value), 2)
    shared_counts <- numeric(ncol(combos))
    shared_lists <- vector("list", ncol(combos))
    combo_names <- character(ncol(combos))
    
    for (i in seq_len(ncol(combos))) {
      a <- combos[1, i]
      b <- combos[2, i]
      combo_name <- paste(a, "&", b)
      shared_feats <- intersect(features_by_value[[a]], features_by_value[[b]])
      
      shared_counts[i] <- length(shared_feats)
      shared_lists[[i]] <- shared_feats
      combo_names[i] <- combo_name
    }
    
    names(shared_counts) <- combo_names
    names(shared_lists) <- combo_names
  }
  
  list(
    value_counts = value_counts,
    features_by_value = features_by_value,
    shared_features_count = shared_counts,
    shared_features_list = shared_lists
  )
}

summary_result <- lapply(multi_list,  function(x) summarize_feature_by_value(x))

# Access outputs
# result <- summary_result[[11]]
# result$value_counts
# result$features_by_value
# result$shared_features_count
# result$shared_features_list

count_signif_allvar <- data.frame()
for(i in 1:length(summary_result)){
  result <- summary_result[[i]]
  mod <- names(summary_result)[i]
  
  signif_list <- unlist(result$features_by_value)
  sum_diet <- length(unique(signif_list[grep("Sago|Pork|Chicken", names(signif_list))]))
  sum_pop<- length(unique(signif_list[grep("DAI|HUL|DPS|PDW", names(signif_list))]))
  sum_age <- length(unique(signif_list[grep("Age", names(signif_list))]))
  sum_sex <- length(unique(signif_list[grep("M", names(signif_list))]))
  sum_bmi <- length(unique(signif_list[grep("bmi", names(signif_list))]))
  
  out <- data.frame(model=mod, Diet=sum_diet, Population=sum_pop, Age=sum_sex, Sex=sum_sex, BMI=sum_bmi)
  
  count_signif_allvar <- rbind(count_signif_allvar, out)
}
count_signif_allvar

bardata <- reshape2::melt(count_signif_allvar, id.vars="model")
bardata$model <- factor(bardata$model, levels(mapvar.long$model))
bardata %>%
  filter(variable %in% c("Diet", "Population")) %>%  # Use filter(), not select(), for row conditions
  group_by(model) %>%
  summarise(total = sum(value), .groups = "drop") -> bartotal
bartotal$variable <- "Diet"
bartotal$variable[grep("population_", bartotal$model)] <- "Population"
bartotal$variable[bartotal$model=="Multivar_diet+age+sex+bmi+population_nonBali"] <- "Diet"
bartotal$model <- factor(bartotal$model, levels(mapvar.long$model))

bardata %>%
  filter(variable %in% c("Age", "Sex","BMI")) %>%  # Use filter(), not select(), for row conditions
  group_by(model) %>%
  summarise(total = sum(value), .groups = "drop") -> barconf
barconf$model <- factor(barconf$model, levels(mapvar.long$model))
barconf$y <- barconf$total+bartotal$total

var_colours <- scales::hue_pal()(5)
names(var_colours) <- c("Diet","Population","Age","Sex","BMI")
  
b <- ggplot(bardata, aes(x=model)) + theme_pubclean()+
  geom_bar(aes(y=value, fill=variable), stat = "identity", position = "stack", width = 0.8) +
  geom_text(data=bartotal, aes(x=model, y=total+30, label=total, colour=variable))+
  geom_text(data=barconf, aes(x=model, y=y+8, label=total), colour="grey50")+
  scale_x_discrete(labels=c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV"))+
  scale_fill_manual(values=var_colours)+
  scale_colour_manual(values=var_colours)+
  scale_y_continuous(limits = c(0,400), expand = c(0,0), breaks = c(5,seq(50,350,by=50)))+
  labs(title = "Count of Microbial Features Associated with Diet and Population",
       subtitle = "Across MaAsLin Model Variants",
       x = NULL, y = "Count") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    #axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "top")
  
b 


# C) Alternative barplot: Blended Variables
base_colours <- c(
  Diet = "#FA8072",        # orange
  Population = "#7A67EE",  # green
  Age = "#76EE00",         # red
  Sex = "#FFD700",         # purple
  BMI = "#40E0D0"          # brown
)

blend_colours <- function(set_names, colour_map, fallback = "#BBBBBB") {
  valid_sets <- intersect(set_names, names(colour_map))
  if (length(valid_sets) == 0) return(fallback)
  rgb_matrix <- sapply(colour_map[valid_sets], function(col) col2rgb(col) / 255)
  if (is.null(dim(rgb_matrix))) return(colour_map[valid_sets])
  rgb_avg <- rowMeans(rgb_matrix)
  rgb(rgb_avg[1], rgb_avg[2], rgb_avg[3])
}

barcount_data <- data.frame()
for(m in mapvar$model){
  #Define list of features per variable
  feature_list <- summary_result[[m]]$features_by_value
  feature_list <- list(Diet=unique(unlist(feature_list[names(feature_list)%in%c("Sago","Pork","Chicken")])),
                       Population=unique(unlist(feature_list[names(feature_list)%in%c("Asmat_DAI","Punan_HUL","Balinese_PDW","Balinese_DPS")])),
                       Age=feature_list$Age,
                       Sex=feature_list$M,
                       BMI=feature_list$bmi)
  
  # Assign each feature to a set combination
  all_features <- unique(unlist(feature_list))
  feature_sets <- map_chr(all_features, function(f) {
    sets <- names(feature_list)[map_lgl(feature_list, ~ f %in% .x)]
    paste(sort(sets), collapse = "_")  # e.g. "Diet_Sex"
  })
  
  # Count how many features per unique overlap group
  overlap_counts <- as.data.frame(table(feature_sets), stringsAsFactors = FALSE)
  colnames(overlap_counts) <- c("overlap", "count")
  overlap_counts$overlap <- as.character(overlap_counts$overlap)
  
  #Assign blended colour for each overlap
  overlap_counts$colour <- sapply(strsplit(overlap_counts$overlap, "_"), function(x) {
    blend_colours(x, base_colours)
  })
  
  overlap_counts <- data.frame(model=m, overlap_counts)
  barcount_data <- rbind(barcount_data, overlap_counts)
}

# Obtain discrete colours 
full_colour <- unique(barcount_data[,c("overlap","colour")])
tmp <- full_colour$colour
names(tmp) <- full_colour$overlap
full_colour <- tmp
full_colour["Population"] <- "#8800aaff"
# Plot
barcount_data$model <- factor(barcount_data$model, levels(mapvar.long$model))
barcount_data$overlap <- factor(barcount_data$overlap, c("Diet","Population","Age","Sex","BMI","BMI_Population","Population_Sex"))
c <- ggplot(barcount_data, aes(x = model, y = count, fill = overlap)) +
  geom_bar(stat = "identity", colour = NA, size = 0.3, width=0.75) +
  scale_fill_manual(values = full_colour) +
  scale_y_continuous(limits = c(0,400), expand = c(0,0), breaks = c(5,seq(50,350,by=50)))+
  scale_x_discrete(labels=c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV"))+
  theme_pubclean() +
  labs(title = "Count of Significant Microbial Features Across MaAsLin Model Variants",
       x = NULL, y = "Count") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    #axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "top")+
  guides(fill = guide_legend(nrow = 1))

c

## compile plot
p3a <- ggarrange(plotlist = list(b+theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 75)),
                          a+theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))), 
          ncol =1,heights = c(3,1))

p3b <- ggarrange(plotlist = list(c+theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 75)),
                          a+theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))), 
          ncol =1,heights = c(3,1))

p3a
p3b

## save plat
#outpath <- "figout/Figure3/"
#for(d in c("svg","png","pdf")){
#  ggsave(plot = p3b, 
#         filename = paste0(outpath,"Fig3D_3.",d), scale = 0.6,
#         device = d,width = 16,height = 8,units = "in",dpi = "print")}
#saveRDS(p3b, file = paste0(outpath,"Fig3D_3.","rds"))


###################
## Age, Sex, BMI ##
###################
multivar.df1 <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/Multivar_diet_popranef/all_results.tsv")
multivar.df2 <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/Multivar_diet+age+sex+bmi_popranef/all_results.tsv")
multivar.df3 <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/Multivar_diet_nonBali_popranef/all_results.tsv")
multivar.df4 <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/Multivar_diet+age+sex+bmi_popranef_nonBali/all_results.tsv")

phenotype_list <- list(base=multivar.df1, adjusted=multivar.df2, base_nonBali=multivar.df3, adjusted_nonBali=multivar.df4)

phenotype_list <- lapply(phenotype_list, function(x) { 
  x <- x[x$qval < 0.05,]
  #x[x$value %in% c("Sago","Pork","Chicken"),"metadata"] <- "Diet"
  return(x)})

phenotype_maglist <- lapply(phenotype_list, function(x) {
  var <- unique(x$value)
  maglist <- lapply(var, function(i) {
    x$feature[x$value == i]
  })
  names(maglist) <- var
  return(maglist)
})

venn_plots <- lapply(phenotype_maglist, function(x) {
  n<-length(x)
  nm <- names(x)
  vdf <- x[nm][1]
  return(vdf)
})

# Remove NULLs (those not plotted due to invalid set size)
# venn_plots <- Filter(Negate(is.null), venn_plots)

##################
## diet volcano ##
##################
df <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/Multivar_diet+age+sex+bmi_popranef/all_results.tsv")
length(unique(df$feature)) -> ct
df.sig <- subset(df, qval < 0.05)
length(unique(df.sig$feature)) -> sig.ct

sig.ct *100 / ct

# append taxonomy
maaslin.taxonomy <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_full_taxonomy.tsv")
maaslin.taxonomy$Species_label <- gsub("_species","_sp.", maaslin.taxonomy$Species_deduplicated) # shorthand taxonomy
maaslin.taxonomy$feature <- maaslin.taxonomy$Genome
df <- merge.data.frame(df, maaslin.taxonomy[,c("feature","original_bin","Genome","Species_label")], by="feature")
colnames(df)[which(colnames(df) == "Species_label")] <- "species"

ggdata <- df
ggdata$association <- ggdata$value
ggdata$association[ggdata$qval>=0.05] <- "none"
ggdata$association[ggdata$association=="M"] <- "Sex:M"
ggdata$association[ggdata$association=="bmi"] <- "BMI"
ggdata$association <- factor(ggdata$association, c("Sago","Pork","Chicken","Age","Sex:M","BMI","none"))

ggdata$abs.coef <- abs(ggdata$coef)
ggdata <- ggdata[order(ggdata$abs.coef,decreasing = T),]
ggdata$label <- ggdata$species
ggdata$label[ggdata$qval>=0.05] <- ""
rownames(ggdata) <- NULL

diet.label <- ggdata[ggdata$label!="",]
diet.label <- merge.data.frame(diet.label, maaslin.taxonomy[,c("feature","Phylum","Family","Species_label")], by="feature")
rownames(diet.label) <- NULL
write.table(diet.label,"output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/all_results_diet_label.tsv",
            quote = F, row.names = F, sep = "\t")

diet.cols <- c("#2E8B57","#FF69B4","#0000FF","#8800aaff","#76EE00","#FFD700","#40E0D0","grey70")
names(diet.cols) <- c("Sago","Pork","Chicken","population","Age","Sex:M","BMI","none")

volc <- ggplot(ggdata, aes(y=-log10(qval), x=coef, fill=association)) + theme_bw() +
  geom_text(aes(label=label), col="grey50", size=3,nudge_y =0.1) +
  annotate("text", x = 0, y = -log10(0.045), label = "qval=0.05", vjust=0, col="red")+
  geom_hline(yintercept = -log10(0.05), col="red", lty=2)+
  geom_point(size=5, pch=21, col="grey30") +
  theme(legend.position = "bottom", legend.box = "horizontal", legend.title = element_blank())+
  guides(fill = guide_legend(nrow = 1)) +
  scale_fill_manual(values = diet.cols) +
  scale_x_continuous(limits = c(-4,4), breaks = sort(c(seq(-4,4, by=1),0)))+
  labs(title = "Diet Correlations (Multi-variate)",
       caption = "q-val is p-val adjusted with Benjamin Hochberg method, q-val<0.05 considered as significant") +
  ylab("Significance -log10(qval)") + xlab("Coefficient")

                               
volc

#outpath <- "figout/diet_maaslin_vocano"
#for(d in c("svg","png","pdf")){
#  ggsave(plot = volc, 
#         filename = paste0(outpath,"diet_volcano.",d), scale = 0.8,
#         device = d,width = 8,height = 8.5,units = "in",dpi = "print")}
#saveRDS(volc, paste0(outpath,"diet_volcano.","rds"))


#######################
## MAASLIN pop ranef ##
#######################
library(lme4)
library(reshape2)
library(dplyr)
library(ggplot2)
library(gridExtra)

# load model 
models <- readRDS("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/Multivar_diet+age+sex+bmi_popranef/fits/models.rds")

partition_fixed_variable=list(Diet=c("Sago","Pork","Chicken"),
                              Age="Age",
                              Sex="GenderM",
                              BMI="bmi")
pop.varex <- data.frame(feature=c(), pop_var=c())
for (i in names(models)){
  test <- models[[i]]
  X <- model.matrix(test)
  beta <- fixef(test)
  eta <- as.matrix(X) %*% beta
  var_fixed <- var(as.numeric(eta))
  
  # partitioned fixed effect
  var_fixed_part <- c()
  for(j in 1:length(partition_fixed_variable)){
    ni <- names(partition_fixed_variable)[j]
    par <- partition_fixed_variable[[j]]
    Xi <- model.matrix(test)[,par]
    bi <- beta[par]
    ei <- as.matrix(Xi) %*% bi
    ivar <- var(as.numeric(ei))
    names(ivar) <- ni
    var_fixed_part <- c(var_fixed_part, ivar)
  }
  
  ranef.var <- data.frame(VarCorr(test))[,c("grp","var1","vcov","sdcor")]
  varex <- ranef.var$vcov
  names(varex) <- ranef.var$grp
  var_random <- varex["population"]
  var_residual <- varex["Residual"]
  
  # proportion of variance
  var_total <- var_fixed + var_random + var_residual
   
  tmp <- data.frame(
    feature = i,
    fixed_var = round(var_fixed * 100 / var_total, 1),
    partial_Diet = round((var_fixed_part["Diet"]/var_fixed) * 100 / var_total, 1),
    partial_Age = round((var_fixed_part["Age"]/var_fixed) * 100 / var_total, 1),
    partial_Sex = round((var_fixed_part["Sex"]/var_fixed) * 100 / var_total, 1),
    partial_BMI = round((var_fixed_part["BMI"]/var_fixed) * 100 / var_total, 1),
    pop_var = round(var_random * 100 / var_total, 1),
    residual_var = round(var_residual * 100 / var_total, 1),
    row.names = NULL
  )
  
  pop.varex <- rbind(pop.varex,tmp)
}

# most abundant in
df.abund <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_comm.tsv",header = T)
df.abund$sampleid <- rownames(df.abund)
df.abund <- melt(df.abund, id.vars = "sampleid",variable.name = "feature",value.name = "CPM")

input.metadata <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_metadata.tsv",header = T)
input.metadata$sampleid <- rownames(input.metadata)
input.metadata <- input.metadata[,c("sampleid","population")]

df.abund <- merge.data.frame(input.metadata,df.abund, by="sampleid")


# append taxonomy
maaslin.taxonomy <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_full_taxonomy.tsv")
maaslin.taxonomy$Species_label <- gsub("_species","_sp.", maaslin.taxonomy$Species_deduplicated) # shorthand taxonomy
maaslin.taxonomy$feature <- maaslin.taxonomy$Genome

df.abund <- merge.data.frame(df.abund, maaslin.taxonomy[,c("feature","original_bin","Genome","Phylum","Species_label")], by="feature")
colnames(df.abund)[which(colnames(df.abund) == "Species_label")] <- "Species"

df.abund %>% group_by(Phylum, Species, feature, population) %>% summarise(popMean=mean(CPM)) -> pop_average

mostabund.pop <- data.frame(feature=c(), topPop=c(), meanCPM=c())
for(i in unique(pop_average$feature)){
  tmp <- pop_average[pop_average$feature==i,]
  max <- max(tmp$popMean)
  n <- grep(max, tmp$popMean)
  
  j <- tmp$population[n]
  sp <- tmp$Species[n]
  ph <- tmp$Phylum[n]
  
  tmp <- data.frame(feature=i, phylum=ph, species=sp, topPop=j, meanCPM=max)
  mostabund.pop <- rbind(mostabund.pop, tmp)
}
mostabund.pop

popranef.df <- merge.data.frame(pop.varex, mostabund.pop, by="feature")
write.table(popranef.df, 
            "output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/Multivar_diet+age+sex+bmi_popranef/all_results_popvar.tsv", 
            quote = F,sep = "\t", row.names = F)

popvar.cutoff <- 35 # in percentage
strong.popranef <- unique(popranef.df[popranef.df$pop_var>=popvar.cutoff,"feature"])
length(strong.popranef)
popranef.out <- popranef.df[popranef.df$feature%in%strong.popranef,]
popranef.out[order(popranef.out$pop_var,decreasing = T),]

######################
## Fig 3E Diet list ##
######################
df <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/Multivar_diet+age+sex+bmi_popranef/all_results.tsv", header=T)
# append taxonomy
maaslin.taxonomy <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_full_taxonomy.tsv")
maaslin.taxonomy$Species_label <- gsub("_species","_sp.", maaslin.taxonomy$Species_deduplicated) # shorthand taxonomy

maaslin.taxonomy$feature <- maaslin.taxonomy$Genome
df <- merge.data.frame(df, maaslin.taxonomy[,c("feature","original_bin","Genome","Phylum","Class","Family","Species_label")], by="feature")
colnames(df)[which(colnames(df) == "Species_label")] <- "species"
colnames(df)[which(colnames(df) == "Phylum")] <- "phylum"
colnames(df)[which(colnames(df) == "Class")] <- "class"
colnames(df)[which(colnames(df) == "Family")] <- "family"

ggdata <- df[,c("feature","value","coef","N.not.0","pval","qval","phylum","class","family","species","Genome")]
ggdata <- ggdata[ggdata$qval<0.05,]
diet.MAGs <- unique(ggdata$feature)
length(diet.MAGs)

# order data
ggdata <- ggdata[order(ggdata$coef,decreasing = F),]
mag.levels <- as.character(unique(ggdata$feature))
ggdata$feature <- factor(ggdata$feature, mag.levels)
phylum.tree.order <- rev(readRDS("output_files/revised_MAG/phylum.tree.order.rds"))
phylum.tree.order <- c(phylum.tree.order,"Methanobacteriota","Thermoplasmatota")
ggdata$phylum <- factor(ggdata$phylum, phylum.tree.order)

# fix label
ggdata$association <- ggdata$value
ggdata$association[ggdata$qval>=0.05] <- "none"
ggdata$association[ggdata$association=="M"] <- "Sex:M"
ggdata$association[ggdata$association=="bmi"] <- "BMI"
ggdata$association <- factor(ggdata$association, c("Sago","Pork","Chicken","Age","Sex:M","BMI","none"))
rownames(ggdata) <- NULL

diet.list.df <- ggdata 

#############################
# abundance heatmap (total) #
#############################
library(reshape2)
df.abund <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_comm.tsv",header = T)
df.abund$sampleid <- rownames(df.abund)
df.abund <- melt(df.abund, id.vars = "sampleid",variable.name = "feature",value.name = "CPM")

input.metadata <- read.delim("output_files/revised_MAG/MAASLIN2/maaslin2_diet+pop_CLR/maaslin_input_metadata.tsv",header = T)
input.metadata$sampleid <- rownames(input.metadata)
input.metadata <- input.metadata[,c("sampleid","population")]

df.abund <- merge.data.frame(input.metadata,df.abund, by="sampleid")
df.abund <- merge.data.frame(df.abund, maaslin.taxonomy[,c("feature","original_bin","Genome","Phylum","Species_label")], by="feature")
colnames(df.abund)[which(colnames(df.abund) == "Species_label")] <- "species"
colnames(df.abund)[which(colnames(df.abund) == "Phylum")] <- "phylum"

n <- length(unique(df.abund$sampleid))

df.abund %>% group_by(feature, phylum) %>% summarise(abundance=sum(CPM)/(n*10^6)) -> total_abund

# order data
mag.levels <- levels(diet.list.df$feature) # make sure to select the correct data frame to filter the features
heat.df.total <- total_abund[total_abund$feature %in% mag.levels,]
heat.df.total$feature <- factor(heat.df.total$feature, mag.levels)
heat.df.total$phylum <- factor(heat.df.total$phylum, phylum.tree.order)


#ggplot(heat.df.total)+ theme_minimal() + 
#  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
#  geom_tile(aes(x="A", y=feature,fill=log(abundance)))+
#  scale_fill_gradient(low = "grey90", high = "black",)+
#  theme(plot.margin = margin(0,0,0,0),
#        panel.spacing = unit(0, "cm"),
#       legend.position = "right",
#      panel.grid = element_blank(),
#     panel.border = element_rect(colour = "black", fill=NA),
#    strip.text.y.left = element_text(angle=0),
#   strip.background = element_rect(fill="white",linewidth = 0),
#  text=element_text(size=12), axis.text=element_text(size=12),
# axis.ticks.x= element_blank(),
#        axis.text.x=element_text(colour = "white"),
#       # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#      axis.text.y = element_text(size=8),
#     axis.title.y= element_blank()) +
#  labs(title="Total Abundance") +
#  scale_x_discrete(expand = c(0,0))+
#  scale_y_discrete(expand = c(0,0))+
#  xlab("")

#############################
# abundance heatmap (bypop) #
#############################
# normalised average
# the normalisation is the percentage of the population average over the sum of the average value in each population by features
df.abund %>% group_by(population,feature, phylum) %>% 
  summarise(n=length(unique(sampleid)),
            popCPM=sum(CPM),
            popMedian=median(CPM),
            popMean=mean(CPM),
            popSD=sd(CPM)) -> pop_abund

df.abund %>% group_by(feature) %>% summarise(totalCPM=sum(CPM)) -> total_abund

heat.df.pop <- merge.data.frame(pop_abund, total_abund, by="feature")

pop_abund %>% group_by(feature) %>% summarise(sum.popMean=sum(popMean)) -> sumpopMean
heat.df.pop <- merge.data.frame(heat.df.pop, sumpopMean, by="feature")

heat.df.pop <- heat.df.pop[order(heat.df.pop$feature),]
heat.df.pop$normalised.abund <- 100*heat.df.pop$popCPM/heat.df.pop$totalCPM
heat.df.pop$normalised.average <- 100*heat.df.pop$popMean/heat.df.pop$sum.popMean

# order data
mag.levels <- levels(diet.list.df$feature) # make sure to select the correct data frame to filter the features
heat.df.pop <- heat.df.pop[heat.df.pop$feature %in% mag.levels,]
heat.df.pop$feature <- factor(heat.df.pop$feature,mag.levels)
heat.df.pop$phylum <- factor(heat.df.pop$phylum, phylum.tree.order)
heat.df.pop$population <- factor(heat.df.pop$population,
                                 c("Asmat_DAI","Punan_HUL","Basap_BRU",
                                   "Balinese_PDW","Balinese_DPS"))

# ggplot(heat.df.pop, aes(x=population, y=feature, fill= normalised.average)) + 
# facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
#theme_minimal() + 
#theme(plot.margin = margin(0,0,0,0),
#        panel.spacing = unit(0, "cm"),
#        legend.position = "right",
#        panel.grid = element_blank(),
#        panel.border = element_rect(colour = "black", fill=NA),
#        text=element_text(size=12), axis.text=element_text(size=12),
#        axis.ticks.x= element_blank(),
#        axis.text.x=element_text(colour = "black", size=8),
#        strip.text.y.left = element_text(angle=0),
#        strip.background = element_rect(fill="white",linewidth = 0),
#        # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#        axis.text.y = element_text(size=8),
#        axis.title.y= element_blank()) +
#  scale_x_discrete(expand = c(0,0))+
#  scale_y_discrete(expand = c(0,0))+
#  geom_tile()+ scale_fill_gradient(low="white",high="black") +
#  labs(title="Normalised Average") +
#  xlab("")

########
# plot #
########
diet.cols <- c("#2E8B57","#FF69B4","#0000FF","#8800aaff","#76EE00","#FFD700","#40E0D0","grey70")
names(diet.cols) <- c("Sago","Pork","Chicken","population","Age","Sex:M","BMI","not significant")

diet.list.df$association <- factor(diet.list.df$association, levels = names(diet.cols))

a <- ggplot(diet.list.df, aes(x=coef,y=feature))+ 
  theme_minimal() + 
  theme(plot.margin = margin(0,0,5,0,unit = "mm"),
        plot.background = element_blank(),
        panel.spacing = unit(0, "cm"),
        panel.grid=element_blank(),
        panel.grid.major.x = element_line(colour="grey",linetype = 3, linewidth = .5),
        panel.grid.minor.x = element_blank(),
        legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size=13),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.y.left = element_text(angle=0, size=10),
        strip.background = element_rect(fill="white",linewidth = 0),
        text=element_text(size=12), axis.text=element_text(size=12),
        #axis.text.y = element_text(angle=0, hjust = 1, size=8),
        axis.text.x=element_text(colour = "black", size=15, color = "grey50"),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y= element_blank())+ 
  guides(fill = guide_legend(title = NA,  nrow = 1, override.aes = list(size = 5)))+
  scale_color_manual(values = diet.cols) +
  scale_x_continuous(breaks = seq(-2,6, by=1),limits = c(-1.5,5.5), expand = c(0.06,0.06))+
  scale_y_discrete(expand = c(0,0.5))+
  xlab("Coefficients") +
  labs(title="MAGs with Significant Diet and Phenotypes Associations (qval<0.05)") +
  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
  geom_point(aes(col=association),size=3) +
  geom_text(aes(label=species), angle=0, nudge_x = 0.1, 
            size=3.5, col="grey40", hjust=0, vjust=0.2)

heat.df.total$feature <- factor(heat.df.total$feature, levels(diet.list.df$feature))
b <- ggplot(heat.df.total)+ 
  theme_minimal() +
  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
  geom_tile(aes(x="Total Abundance", y=feature,fill=log(abundance)))+
  scale_fill_gradient(low = "grey90", high = "black",)+
  theme(plot.margin = margin(0,0,0,0,unit = "mm"),
        plot.background = element_blank(),
        panel.spacing = unit(0, "cm"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        #strip.text.y.left = element_text(angle=0),
        #strip.background = element_rect(fill="white",linewidth = 0),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(size=10, colour = "black",angle=60, hjust = 1),
        axis.ticks.x= element_line(colour = NA),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        # axis.text.y = element_text(size=8),
        axis.title.y= element_blank()) +
  labs(title="") +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0.5))+
  xlab("")

heat.df.pop$feature <- factor(heat.df.pop$feature, levels(diet.list.df$feature))
c <- ggplot(heat.df.pop, aes(x=population, y=feature, fill=normalised.average)) + 
  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
  theme_minimal() + 
  theme(plot.margin = margin(0,0,2,0,unit = "mm"),
        plot.background = element_blank(),
        panel.spacing = unit(0, "cm"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.grid = element_blank(),
        #strip.text.y.left = element_text(angle=0),
        #strip.background = element_rect(fill="white",linewidth = 0),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(size=10, colour = "black",angle=60, hjust = 1),
        axis.ticks.x= element_line(colour = NA),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        # axis.text.y = element_text(size=8),
        axis.title.y= element_blank()) +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0.5))+
  geom_tile()+ scale_fill_gradient(low="grey95",high="black") +
  labs(title="") +
  xlab("")

p <- ggarrange(plotlist = list(a,b,c), nrow = 1,widths = c(22,1,4))
p

# draw plot
outpath <- "figout/Figure3/"
for(d in c("svg","png","pdf")){
  ggsave(plot = p, 
         filename = paste0(outpath,"Fig3E.",d), scale = 0.8,
         device = d,width = 9,height = 10,units = "in",dpi = "print")}
saveRDS(p, file = paste0(outpath,"Fig3E.","rds"))

# print duplicated
diet.list.df[duplicated(diet.list.df$feature),] # remember to post-edit if length != 0

##################
### POP RANEF ####
##################
# create figure data frame
ggdata <- popranef.df
ggdata$topPop <- factor(ggdata$topPop, c("Asmat_DAI","Basap_BRU","Punan_HUL","Balinese_DPS", "Balinese_PDW","none"))
ggdata <- ggdata[ggdata$pop_var!=0,]
ggdata$label <- ggdata$feature

summary(ggdata$pop_var)
quantile(ggdata$pop_var)
hist(ggdata$pop_var, breaks=30)
quantile(ggdata$pop_var)

ggdata %>% summarise(n=length(pop_var), 
                     average=mean(pop_var),
                     sd=sd(pop_var),
                     se=sd(pop_var)/length(pop_var),
                     min=min(pop_var),
                     Q1=quantile(pop_var)[2],
                     med=median(pop_var),
                     Q3=quantile(pop_var)[4],
                     max=max(pop_var))

ggdata %>% 
  filter(pop_var > 35) %>% 
  group_by(phylum) %>%
  summarise(n=length(pop_var), 
            average=mean(pop_var),
            sd=sd(pop_var),
            se=sd(pop_var)/length(pop_var),
            min=min(pop_var),
            Q1=quantile(pop_var)[2],
            med=median(pop_var),
            Q3=quantile(pop_var)[4],
            max=max(pop_var))



# Filter MAGs 
ggdata <- ggdata[order(ggdata$pop_var, decreasing = T),]
best.pop <- ggdata[ggdata$pop_var>35,] # varex > 35%
# best.pop <- ggdata 

# set phylum colours
phylum.cols <- read.delim("input_files/phylum_colours_bar.tsv", header=T)
tmp <- phylum.cols$cols
names(tmp) <- phylum.cols$phylum
phylum.cols <- tmp

# set population colours
popcols <- readRDS("input_files/pop_mycols.rds")
popcols <- popcols[names(popcols) %in% unique(input.metadata$population)]

table(best.pop$phylum)

# order data
best.pop <- best.pop[order(best.pop$pop_var,decreasing = T),]
best.pop$feature <- factor(best.pop$feature, level=unique(best.pop$feature))
mag.levels <- unique(best.pop$feature)
best.pop$feature <- factor(best.pop$feature, rev(mag.levels))
best.pop <- best.pop[order(rev(mag.levels)),]

phylum.tree.order <- rev(readRDS("output_files/revised_MAG/phylum.tree.order.rds"))
phylum.tree.order <- c(phylum.tree.order,"Methanobacteriota","Thermoplasmatota")
best.pop$phylum <- factor(best.pop$phylum, phylum.tree.order)

rownames(best.pop) <- NULL

# add diet associations
diet.assoc <- df[df$qval<0.05, c("feature","value")]
best.pop$diet.assoc <- "population"
for(i in 1:nrow(best.pop)){
  x <- as.character(best.pop$feature)[i]
  y <- x %in% diet.assoc$feature
  z <- diet.assoc$value[diet.assoc$feature==x]
  if(isTRUE(y)){
    best.pop[i,"diet.assoc"] <- y
  }
}

diet.cols <- c("#2E8B57","#FF69B4","#0000FF","#8800aaff","#76EE00","#FFD700","#40E0D0","grey70")
names(diet.cols) <- c("Sago","Pork","Chicken","population","Age","Sex:M","BMI","not significant")

#ggplot(best.pop, aes(x=pop_var,y=feature))+ theme_classic() +
#  theme(plot.margin = margin(0,0,0,0,unit = "mm"),
#        panel.spacing = unit(0, "cm"),
#        legend.position = "right",
#        panel.background = element_rect(color = "black",fill=NA),
#        strip.text.y.left = element_text(angle=0),
#        strip.background = element_rect(fill="white",linewidth = 0),
#        text=element_text(size=12), axis.text=element_text(size=12),
#        axis.text.y = element_text(angle=0, hjust = 1, size=8),
#        axis.text.x=element_text(colour = "black", size=8),
#        # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#        axis.title.y= element_blank(),
#        panel.grid.major.x = element_line(colour="grey",linetype = 2)) + 
#  scale_color_manual(values = diet.cols) +
#  scale_x_continuous(breaks = seq(0,100, by=10),limits = c(35,105), expand = c(0,0),labels = paste0(seq(0,100, by=10),"%"))+
#  scale_y_discrete(expand = c(0,0.5))+
#  xlab("Variance Explained") +
#  labs(title="MAGs with Strongest Population Effect (varience explained >35%)")+
#  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
#  geom_point(aes(col=diet.assoc),size=3) +
#  geom_text(aes(label=species), angle=0, nudge_x = 0.5, 
#            size=3, col="grey40", hjust=0, vjust=0.2)

###################
# Total Abundance #
###################
n <- length(unique(df.abund$sampleid))
df.abund %>% group_by(feature, phylum) %>% summarise(abundance=sum(CPM)/(n*10^6)) -> total_abund

# order data
mag.levels <- as.character(unique(best.pop$feature))
heat.df.total <- total_abund[total_abund$feature %in% as.character(mag.levels),]
heat.df.total$feature <- factor(heat.df.total$feature, mag.levels)
heat.df.total$phylum <- factor(heat.df.total$phylum, phylum.tree.order)

# ggplot(heat.df.total)+ theme_classic() +
#  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
#  geom_tile(aes(x="%Abundance (log10)", y=feature,fill=log(abundance)))+
#  scale_fill_gradient(low = "grey90", high = "black",)+
#  theme(plot.margin = margin(0,0,0,0,unit = "mm"),
#        panel.spacing = unit(0, "cm"),
#        legend.position = "right",
#        panel.grid = element_blank(),
#        panel.border = element_rect(colour = "black", fill=NA),
#        strip.text.y.left = element_text(angle=0),
#        strip.background = element_rect(fill="white",linewidth = 0),
#        text=element_text(size=12), axis.text=element_text(size=12),
#        axis.ticks.x= element_blank(),
#        axis.text.x=element_text(colour = "black", size=8),
#        # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#        axis.text.y = element_text(size=8),
#        axis.ticks.y= element_blank(),
#        axis.title.y= element_blank()) +
#  labs(title="Total Abundance") +
#  scale_x_discrete(expand = c(0,0))+
#  scale_y_discrete(expand = c(0,0.5))+
#   xlab("")

########################
# population abundance #
########################
# normalised average
df.abund %>% group_by(population,feature, phylum) %>% 
  summarise(popCPM=sum(CPM),
            popMedian=median(CPM),
            popMean=mean(CPM),
            popSD=sd(CPM)) -> pop_abund

df.abund %>% group_by(feature) %>% summarise(totalCPM=sum(CPM)) -> total_abund

heat.df.pop <- merge.data.frame(pop_abund, total_abund, by="feature")

pop_abund %>% group_by(feature) %>% summarise(sum.popMean=sum(popMean)) -> sumpopMean
heat.df.pop <- merge.data.frame(heat.df.pop, sumpopMean, by="feature")

heat.df.pop <- heat.df.pop[order(heat.df.pop$feature),]

heat.df.pop$normalised.abund <- 100*heat.df.pop$popCPM/heat.df.pop$totalCPM
heat.df.pop$normalised.average <- 100*heat.df.pop$popMean/heat.df.pop$sum.popMean

# order data
mag.levels <- as.character(unique(best.pop$feature))
heat.df.pop <- heat.df.pop[heat.df.pop$feature %in% as.character(mag.levels),]
heat.df.pop$feature <- factor(heat.df.pop$feature,mag.levels)
heat.df.pop$phylum <- factor(heat.df.pop$phylum, phylum.tree.order)
heat.df.pop$population <- factor(heat.df.pop$population,
                                 c("Asmat_DAI","Punan_HUL","Basap_BRU",
                                   "Balinese_PDW","Balinese_DPS"))

# ggplot(heat.df.pop, aes(x=population, y=feature, fill=normalised.average)) + 
#  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
#  theme_minimal() + 
#  theme(plot.margin = margin(0,0,0,0,unit = "mm"),
#        panel.spacing = unit(0, "cm"),
#        legend.position = "right",
#        panel.grid = element_blank(),
#        panel.border = element_rect(colour = "black", fill=NA),
#        strip.text.y.left = element_text(angle=0),
#        strip.background = element_rect(fill="white",linewidth = 0),
#        text=element_text(size=12), axis.text=element_text(size=12),
#        axis.ticks.x= element_blank(),
#        axis.text.x=element_text(colour = "black", size=8),
#        # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#        axis.text.y = element_text(size=8),
#        axis.ticks.y= element_blank(),
#        axis.title.y= element_blank()) +
#  scale_x_discrete(expand = c(0,0))+
#  scale_y_discrete(expand = c(0,0.5))+
#  geom_tile()+ scale_fill_gradient(low="white",high="black") +
#  labs(title="Normalised Average") +
#  xlab("")

###################
# PLOT (Popranef) #
###################
a <- ggplot(best.pop, aes(x=pop_var,y=feature,unit = "px"))+ 
  theme_minimal() + 
  theme(plot.margin = margin(0,0,14.5,0, unit = "mm"),
        plot.background = element_blank(),
        panel.spacing = unit(0, "cm"),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour="grey",linetype = 3, linewidth = .5),
        panel.grid.minor.x = element_line(colour="grey",linetype = 3, linewidth = .5),
        legend.position = "none",
        panel.background = element_rect(fill="white"),
        panel.border = element_rect(colour = "black",fill = NA),
        strip.text.y.left = element_text(angle=0,size = 10),
        strip.background = element_rect(fill="white",linewidth = 0),
        text=element_text(size=12), axis.text=element_text(size=12),
        #axis.text.y = element_text(angle=0, hjust = 1, size=8),
        axis.text.x=element_text(colour = "grey50", size=15),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y= element_blank())+
  scale_color_manual(values = diet.cols) +
  scale_x_continuous(breaks = seq(0,100, by=20),limits = c(35,110), 
                     expand = c(0.06,0.06),labels = paste0(seq(0,100, by=20),"%"))+
  scale_y_discrete(expand = c(0,0.5))+
  xlab("Variance Explained") +
  labs(title="MAGs with Population Effect (variance explained >35%)") +
  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
  geom_point(aes(col=diet.assoc),size=3)+
  geom_text(aes(label=species), angle=0, nudge_x = 2, 
            size=3.5, col="grey40", hjust=0, vjust=0.5)

heat.df.total$feature <- factor(heat.df.total$feature, levels(best.pop$feature))
b <- ggplot(heat.df.total)+ 
  theme_minimal() + 
  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
  geom_tile(aes(x="Total Abundance", y=feature,fill=log(abundance)))+
  scale_fill_gradient(low = "grey90", high = "black",)+
  theme(plot.margin = margin(0,0,0,0,unit = "mm"),
        plot.background = element_blank(),
        panel.spacing = unit(0, "cm"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        # strip.text.y.left = element_text(angle=0),
        # strip.background = element_rect(fill="white",linewidth = 0),
        strip.text = element_blank(),
        strip.background = element_blank(),
        text=element_text(size=12), axis.text=element_text(size=12),
        axis.ticks.x= element_blank(),
        axis.text.x=element_text(colour = "black", size=8, angle = 60, hjust = 1),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        #axis.text.y = element_text(size=8), axis.ticks.y= element_blank(),
        axis.title.y= element_blank()) +
  labs(title="") +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0.5))+
  xlab("")

heat.df.pop$feature <- factor(heat.df.pop$feature, levels(best.pop$feature))
c <- ggplot(heat.df.pop, aes(x=population, y=feature, fill=normalised.average)) + 
  facet_grid(phylum~.,switch = "y",space = "free_y", scales = "free_y") +
  theme_minimal() + 
  theme(plot.margin = margin(0,0,1.8,0,unit = "mm"),
        plot.background = element_blank(),
        panel.spacing = unit(0, "cm"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        # strip.text.y.left = element_text(angle=0),
        # strip.background = element_rect(fill="white",linewidth = 0),
        strip.text = element_blank(),
        strip.background = element_blank(),
        text=element_text(size=12), axis.text=element_text(size=12),
        axis.ticks.x= element_blank(),
        axis.text.x=element_text(colour = "black", size=8,angle = 60, hjust = 1),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        # axis.text.y = element_text(size=8),axis.ticks.y= element_blank(),
        axis.title.y= element_blank()) +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0.5))+
  geom_tile()+ scale_fill_gradient(low="grey95",high="black") +
  labs(title="") +
  xlab("")

q <- ggarrange(plotlist = list(a,b,c), nrow = 1,widths = c(20,1,4))
q

# draw plot
outpath <- "figout/Figure3/"
for(d in c("svg","png","pdf")){
  ggsave(plot = q, 
         filename = paste0(outpath,"Fig3F.",d), scale = 0.8,
         device = d,width = 9,height = 12,units = "in",dpi = "print")}
saveRDS(q, file = paste0(outpath,"Fig3F.","rds"))


############################################
# Average pop varex in diet-associated MAGs #
############################################
diet.popvarex <- popranef.df
diet.popvarex$significant <- "not_significant"
diet.popvarex$significant[diet.popvarex$feature%in%diet.MAGs] <- "significant"

diet.popvarex <- melt(diet.popvarex[,c("feature","significant","fixed_var","pop_var","residual_var")], id.vars = c("feature","significant"), variable.name = "variance")

ggplot(diet.popvarex, aes(x=significant, y=value, colour=variance))+ 
  theme_classic() +
  geom_point(size=2, position = position_jitterdodge(jitter.width = 0.2,dodge.width = 0.6,jitter.height = 0)) +
  geom_boxplot(outlier.shape = NA, alpha=0.5, fill="white", position = position_dodge(width = 0.6))

psych::describe(popranef.df[popranef.df$feature%in%diet.MAGs, "pop_var"])
psych::describe(popranef.df[!popranef.df$feature%in%diet.MAGs, "pop_var"])

wilcox.test(diet.popvarex[diet.popvarex$variance=="pop_var","value"] ~ diet.popvarex[diet.popvarex$variance=="pop_var","significant"])

