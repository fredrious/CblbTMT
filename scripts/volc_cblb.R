

volc_cblb <- function(result, NominalCutoff=1, AdjustedCutoff=0.05, LabellingCutoff=0.05, FCCutoff=1, savefile=savefile, label=TRUE) {
  
  library(ggrepel)
  
  # limma_top <- setDT(read.csv(file="results/cblb_norm/revised/topTable_limma.csv", header=TRUE, sep=","))
  # limma_top$Comparison <- gsub("group", "", limma_top$Comparison)

  # AdjustedCutoff=0.05
  # FCCutoff=1
  # LabellingCutoff=0.05
  # NominalCutoff=1
  
  setDT(result)
  result$Significance <- "NS"
  result$Significance[(abs(result$logFC) > FCCutoff)] <- "FC"
  result$Significance[(result$adj.P.Val < AdjustedCutoff)] <- "P.adj"
  result$Significance[(result$adj.P.Val < AdjustedCutoff) & (abs(result$logFC) > FCCutoff)] <- "FC_P.adj"
  table(result$Significance)
  result$Significance <- factor(result$Significance, levels=c("NS", "FC", "P.adj", "FC_P.adj"))
  result$Label <- as.factor(result$Comparison)
  
  
  
  pdf(savefile)
  for (ii in 1:uniqueN(result$Label)) {
    
    sub <- result[Label == unique(result$Label)[ii]]
    mainn <- as.character(unique(sub$Label))
    
    # limma_top_sub <- limma_top[Comparison %in% main, list(Protein, rank=as.numeric(seq(X))Z][1:50,] #top 50 proteins by limma per comparison
    # merge(sub,limma_top_sub, all.x=TRUE)
    
    plot2 <- ggplot(sub, aes(x=logFC, y=-log10(P.Value))) +
      
      #Add points:
      #	Colour based on factors set a few lines up
      #	'alpha' provides gradual shading of colour
      #	Set size of points
      geom_point(aes(color=factor(Significance)), alpha=1/2, size=0.8) +
      
      #Choose which colours to use; otherwise, ggplot2 choose automatically 
      # (order depends on how factors are ordered in result$Significance)
      scale_color_manual(values=c(NS="grey30", FC="forestgreen", P.adj="royalblue", FC_P.adj="red2"), 
                         labels=c(NS="NS", 
                                  FC=paste("LogFC>|", FCCutoff, "|", sep=""), 
                                  P.adj=paste("P.adj<", AdjustedCutoff, sep=""), 
                                  FC_P.adj=paste("P.adj<", AdjustedCutoff, " & LogFC>|", FCCutoff, "|", sep=""))) +
      
      #Set the size of the plotting window
      theme_bw(base_size=24) +
      # facet_grid(first~second) +
      # facet_wrap(~Label) +
      
      
      #Modify various aspects of the plot text and legend
      theme(legend.background=element_rect(),
            plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
            
            panel.grid.major=element_blank(),	#Remove gridlines
            panel.grid.minor=element_blank(),	#Remove gridlines
            
            axis.text.x=element_text(angle=0, size=12, vjust=1),
            axis.text.y=element_text(angle=0, size=12, vjust=1),
            axis.title=element_text(size=12),
            strip.text.x = element_text(size=10, margin = margin(1,0,1,0)),
            
            #Legend
            legend.position="top",			#Moves the legend to the top of the plot
            legend.key=element_blank(),		#removes the border
            legend.key.size=unit(0.5, "cm"),	#Sets overall area/size of the legend
            legend.text=element_text(size=8),	#Text size
            title=element_text(size=8),		#Title text size
            legend.title=element_blank()) +		#Remove the title
      
      #Change the size of the icons/symbols in the legend
      guides(colour = guide_legend(override.aes=list(size=2.5))) +
      
      #Set x- and y-axes labels
      xlab(bquote(~Log[2]~ "fold change")) +
      ylab(bquote(~-Log[10]~italic(P))) +
      
      #Set the axis limits
      #xlim(-6.5, 6.5) +
      #ylim(0, 100) +
      
      #Add a vertical line for fold change cut-offs
      geom_vline(xintercept=c(-FCCutoff, FCCutoff), linetype="longdash", colour="black", size=0.4) +
      geom_vline(xintercept=0, linetype="dotted", colour="grey", size=0.4) +
      
      #Add a horizontal line for P-value cut-off
      geom_hline(yintercept=-log10(AdjustedCutoff), linetype="longdash", colour="black", size=0.4) +
      #Set title
      ggtitle(mainn) 
      
      if (label) {
      plot2 <- plot2 +
      #Tidy the text labels for a subset of genes
      geom_text_repel(data=subset(sub[Gene=="Cblb"], P.Value < LabellingCutoff & abs(logFC) > FCCutoff),
                # aes(label=rownames(subset(result, adj.pvalue < LabellingCutoff & abs(log2FC) > FCCutoff))),
                aes(label=paste(Protein,Gene, sep = "\n")),
                size=2.3,
                #segment.color="black", #This and the next parameter spread out the labels and join them to their points by a line
                #segment.size=0.01,
                # check_overlap=TRUE,
                vjust=0)
      
      }
    
    
    
    
    nm <- gsub("\\.","",as.character(unique(sub$Label)))
    # ggsave(paste0(nm,".pdf"), width = 20, height = 20, units = "cm")
    print(plot2)
  }
  dev.off()
}
