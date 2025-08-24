mypop <- read.delim("input_files/mypop_wSPMP.tsv")
mypop$urbanisation <- factor(mypop$urbanisation,c("remote","rural","urban"))
mypop <- mypop[order(mypop$urbanisation),]
mypop$population <- factor(mypop$population, levels = readRDS("input_files/population_order.rds"))
scount <- data.frame(table(mypop$population), stringsAsFactors = F)
colnames(scount) <- c("population","n")
scount$group <- paste0(scount$population,"\n(n=", scount$n,")")
mypop <- merge.data.frame(mypop, scount, by="population")
rm(scount)
saveRDS(mypop,"input_files/mypop_wSPMP.rds")
mypop <- readRDS("input_files/mypop_wSPMP.rds")
pop.ord <- levels(mypop$population)

# colours
mycols <- read.delim("input_files/MAGs_pop_mycols.tsv", header=T)
mycols <- mycols[mycols$group %in% unique(mypop$population),]
col.ord <- mycols$group
mycols <- mycols$mycols
names(mycols) <- col.ord
mycols <- mycols[pop.ord]

# shapes
myshape <- c(21,22,24)
names(myshape) <- c("remote","rural","urban")