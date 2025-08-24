# compare families
sampleid <- c(readRDS("input_files/sampleid_n116.rds"), readRDS("input_files/sampleid_spore.rds"))
mypop <- read.delim("input_files/mypop_wSPMP.tsv")
mypop <- mypop[mypop$sampleid%in% sampleid,]
myclust <- stats::cutree(readRDS("figout/additional_figs/ward_clusters_SIN.rds"), 3)
mypop$clust_k3 <- myclust[mypop$sampleid]
mypop$clust_k3[mypop$clust_k3 == 1] <- "Cluster"
mypop$clust_k3[mypop$clust_k3 == 3] <- "Cluster2" # flipped to match ls transition order
mypop$clust_k3[mypop$clust_k3 == 2] <- "Cluster3"
table(mypop$population, mypop$clust_k3)

full.taxonomy <- read.delim("input_files/taxonomy/indomee_full_taxonomy.tsv")
colnames(full.taxonomy)[3:10] <- tolower(colnames(full.taxonomy)[3:10])

relab <- read.delim("input_files/abundance_tables_raw/original/bwa_counts_total_filtered_wMetadata_HybridAssembly_Comp50Cont5_draft2_strongFilter_relAbund_revised.tsv")
RA <- data.frame(relab[,sampleid], row.names=relab$Genome)
colSums(RA)
RA <- round(RA * 10^6,0)
RA$Genome <- rownames(RA)
RA <- reshape2::melt(RA, id.vars="Genome", variable.name="sampleid", value.name="CPM")
RA <- merge.data.frame(RA, full.taxonomy[,c("Genome","family")])
RA <- merge.data.frame(RA, mypop[,c("sampleid","population","urbanisation","clust_k3")])
summary(RA$CPM[RA$CPM!=0])

library(dplyr)
family_abund <- RA %>% group_by(sampleid, population, urbanisation, clust_k3, family) %>% summarise(CPM=sum(CPM)) 
head(family_abund)
summary(family_abund$CPM[family_abund$CPM!=0])

family.kruskal.out <- do.call(rbind,lapply(unique(family_abund$family), function(x) {
  tmp <- subset(family_abund, family == x)
  out <- kruskal.test(tmp$CPM ~ tmp$clust_k3)
  out <- data.frame(family=x, test=out$method, kruskal_chisq=round(as.numeric(out$statistic),2),
                    N = length(tmp$CPM), k=length(unique(tmp$clust_k3)),
                    pval=out$p.value, qval=p.adjust(out$p.value, n = length(unique(family_abund$family))))
  out$effect_size <- (out$kruskal_chisq - out$k + 1) / (out$N - out$k) # eta square
  out <- out[, c(1,2,3,4,5,8,6,7)]
  out <- cbind(out, data.frame(meanCPM=round(mean(tmp$CPM),0), totalCPM=sum(tmp$CPM), 
                               stderrCPM=round((sd(tmp$CPM)/sqrt(length(tmp$CPM))),0)))
  return(out)
}))

family.kruskal.out <- merge.data.frame(unique(full.taxonomy[,c("phylum","class","order","family")]),
                                       family.kruskal.out)
family.kruskal.out <- family.kruskal.out[order(family.kruskal.out$effect_size, decreasing = T),]

head(family.kruskal.out)
sum(family.kruskal.out$pval<0.05)
sum(family.kruskal.out$qval<0.05) # adjusted
round(100* sum(family.kruskal.out$qval<0.05) / nrow(family.kruskal.out)) # % significant (adjusted)

clustMean_RA <- RA %>% group_by(family, ifelse(clust_k3=="Cluster3","Cluster3","Cluster1+2")) %>% 
  summarise(n=length(CPM),meanCPM=round(mean(CPM),0), logMean=log10(mean(CPM)))
colnames(clustMean_RA)[2] <- "group"
head(clustMean_RA)

tmp <- reshape2::recast(family ~ group, data=clustMean_RA[,c("family","group","meanCPM")])
family.kruskal.out <- merge.data.frame(family.kruskal.out, tmp)
family.kruskal.out$logRatio <-  log10(family.kruskal.out$"Cluster3") / log10(family.kruskal.out$"Cluster1+2")
family.kruskal.out$log2_fold_change <- log2((family.kruskal.out$"Cluster3" + 1) / (family.kruskal.out$"Cluster1+2" + 1))
family.kruskal.out$fold_change <- sign(family.kruskal.out$log2_fold_change) * (2^abs(family.kruskal.out$log2_fold_change))

family.kruskal.out$occur <- "Both"
family.kruskal.out[family.kruskal.out$Cluster3==0,"occur"] <- "Cluster1+2"
family.kruskal.out[family.kruskal.out$"Cluster1+2"==0,"occur"] <- "Cluster3"
colnames(family.kruskal.out)[c(15,16)] <- paste0("meanCPM_",colnames(family.kruskal.out)[c(15,16)])

head(family.kruskal.out)
write.table(family.kruskal.out, file = "figout/additional_figs/family_vsCluster_kruskal.tsv", quote = F,sep = "\t",row.names = F)

# Bifidobacteriaceae
family.kruskal.out[family.kruskal.out$family=="Bifidobacteriaceae",]
