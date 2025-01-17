# require(ggplot2)
# require(ggpubr)
# require(gridExtra)
# require(vegan)
axis_labels <- function(cap){
  eig <- cap$CA$eig
  var <- eig/sum(eig)
  axis_labs <- paste("PCo",1:length(names(var)), " (", round(var*100,1),"%)", sep="")
  axis_labs
}

pretty_PCA <- function(cap, 
                       axes, 
                       metadata, 
                       colour.by, 
                       colour.scheme,
                       colour.order,
                       shape.by,
                       shape.scheme,
                       shape.order,
                       output.path=NULL){
    # select PC
    PC <- axes
    
    # prepare plot data for collation
    xycoord <- summary(cap)
    xycoord <- data.frame(xycoord$sites[,PC])
    xycoord$sampleid <- rownames(xycoord)
    rownames(xycoord) <- NULL
    
    colnames(xycoord) <- c("Xaxis","Yaxis","sampleid")
    
    # axes limits
    xlims <- c(min(xycoord$Xaxis), max(xycoord$Xaxis)) * 1.1
    ylims <- c(min(xycoord$Yaxis), max(xycoord$Yaxis)) * 1.1
    
    # set environmental variables
    ggdata <- merge.data.frame(xycoord, metadata, by="sampleid")
    ggdata <- ggdata[,c(colnames(xycoord), colour.by, shape.by)]
    colnames(ggdata) <- c(colnames(xycoord),"my_colour","my_shape")
    
    # set sample order
    pop.ord <- colour.order
    
    ggdata$my_colour <- factor(ggdata$my_colour, pop.ord)
    ggdata <- ggdata[order(ggdata$my_colour),]
    ggdata$sampleid <- factor(ggdata$sampleid, unique(ggdata$sampleid))
    rownames(ggdata) <-  NULL
    
    # calculate eigen values and set axis title
    eig <- cap$CA$eig
    var <- eig/sum(eig)
    axis_labs <- paste("PCo",1:length(names(var)), " (", round(var*100,1),"%)", sep="")
    xtitle <- axis_labs[PC[1]]
    ytitle <- axis_labs[PC[2]]
    
    # order colour scheme
    mycols <- colour.scheme[pop.ord]
    
    # order shape scheme
    shape.scheme <- shape.scheme[shape.order]
    
    # make base plot
    a <- ggplot() +
      geom_vline(xintercept = 0, lty="dashed")+
      geom_hline(yintercept=0, lty="dashed")
    
    # draw the rest of the plot
    a <- a + 
      theme_bw(base_size = 12) +
      geom_point(data=ggdata, aes(x=Xaxis, y=Yaxis, fill=my_colour, shape=my_shape, size=my_shape), col="black") + 
      scale_fill_manual(values=mycols) +
      scale_shape_manual(values=shape.scheme)+
      scale_size_manual(values=c(3.5,3.2,2.7)) +
      scale_x_continuous(expand = c(0,0),
                         limits = xlims,
                         breaks = seq(-10,10,1),
                         position = "top")+
      scale_y_continuous(expand = c(0,0),
                         limits = ylims,
                         breaks = seq(-10,10,1),
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
    ggdata$my_colour <- factor(ggdata$my_colour, rev(pop.ord))
    ggdata <- ggdata[order(ggdata$my_colour),]
    b <- ggplot(ggdata,
                aes(x=my_colour, y=Xaxis, fill=my_colour))+
      theme_bw(base_size = 12) +
      geom_hline(yintercept = 0, lty="dashed") +
      geom_boxplot(width=0.8, outlier.size = 2,
                   outlier.shape=21, outlier.fill = "white") +
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
                         breaks = seq(-10,10,1),
                         position = "left") +
      xlab("sample distance") + ylab(" ")+
      coord_flip()
    
    # PCo Group Boxplot Y axis
    ggdata$my_colour <- factor(ggdata$my_colour, pop.ord)
    c <- ggplot(ggdata,
                aes(x=my_colour, y=Yaxis, fill=my_colour)) +
      theme_bw(base_size = 12) +
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
                         breaks = seq(-10,10,1),
                         position = "right")+
      scale_x_discrete(position="top")+
      xlab("sample distance") + ylab(" ")
    
    # fix labels and extract legend
    pop_count <- data.frame(table(metadata[,colour.by]))
    colnames(pop_count)  <- c(colour.by,"sample")
    pop_count[,colour.by] <- factor(pop_count[,colour.by], pop.ord)
    pop_count
    
    poplab <-paste0(pop_count[,colour.by]," (n=",pop_count$sample,")")
    names(poplab) <- pop_count[,colour.by]
    
    ggdata$label <- c()
    for(pop in names(poplab)){
      ggdata[ggdata$my_colour == pop, "label"] <- poplab[pop]
    }
    ggdata <- ggdata[order(ggdata$my_colour),]
    ggdata$label <- factor(ggdata$label, unique(ggdata$label))
    
    legend.col <-  mycols[as.character(unique(ggdata$my_colour))]
    names(legend.col) <- poplab[names(legend.col)]
    
    p <- ggplot(ggdata, aes(x=label, y=Xaxis, fill=label))+ theme_bw(base_size = 10) +
      geom_boxplot(width=0.8, outlier.size = 2, outlier.shape=21, outlier.fill = "white") +
      scale_fill_manual(values=legend.col) +
      theme(legend.text = element_text(size=10),
            legend.title = element_blank(),
            legend.key.size = unit(0.3,"cm"),
            plot.background = element_rect(colour = "white"))
    
    p <- ggplot(ggdata)+ theme_void(base_size = 5) +
      geom_point(aes(x=label, y=Xaxis, shape=my_shape))+
      geom_boxplot(aes(x=label, y=Xaxis, fill=label))+    
      scale_fill_manual(values=legend.col) +
      scale_shape_manual(values=shape.scheme)+
      theme(legend.text = element_text(size=5),
            legend.title = element_blank(),
            legend.key.size = unit(0.3,"cm"),
            legend.box.margin = margin(0,0,0,0,"cm"))
    
    legend <- cowplot::get_legend(p)
    legend <- as_ggplot(legend) + theme(plot.margin = margin(-1.5,1,0,0,"cm"))
    
    
    # set plot layout
    lay <- rbind(c(1,1,1,1,2,2),
                 c(1,1,1,1,2,2),
                 c(1,1,1,1,2,2),
                 c(1,1,1,1,2,2),
                 c(3,3,3,3,4,4),
                 c(3,3,3,3,4,4))
    
    
    if(is.null(output.path)){
      # draw plot
      grid.arrange(a,c,b,legend, layout_matrix=lay)
    }else{
    #save.plot
    filename=paste(paste0("PCo",c(1,2)), collapse = "_")
    ggsave(plot=grid.arrange(a,c,b,legend, layout_matrix=lay),
           filename = paste0(output.path,filename,".png"), scale = 1.8,
           device = "png",dpi = "print" ,width = 10, height = 10,units = "cm")
    
    ggsave(plot=grid.arrange(a,c,b,legend, layout_matrix=lay),
           filename = paste0(output.path,filename,".svg"), scale = 1.8,
           device = "svg",dpi = "print" ,width = 10, height = 10,units = "cm")
    
    print(paste0("Output files saved to: ",output.path))
    }
}
       
#usage
# pretty_PCA(cap = ordi, metadata = mypop, 
  #         colour.by = "population", colour.scheme = mycols,axes = c(1,2),
   #        path.to.save=mypath, filename = myfile)
