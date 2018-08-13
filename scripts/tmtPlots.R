

################################
## function to plot peptides
ppt.plot <- function(ddt) {
  tmp <- data.table(Peptide = as.character(unique(ddt$Peptide)), PeptideN=NA)
  tmp[, PeptideN := sapply(Peptide, function(x) .spltpep(x,n=9))]
  ddt <- ddt[tmp, on="Peptide"]
  
  p.tile <- ggplot(data=ddt, aes(x=Run, y=PeptideN, fill=Abundance)) + 
    geom_tile(color="white", size=0.05) +
    scale_fill_distiller(palette = "Spectral",na.value = "grey70") +
    facet_grid(PeptideN + Charge ~ Mixture+Channel, scales = "free", space = "free", switch="y") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.text.y=element_blank()) +
    ggtitle(paste0("Prot: ", unique(ddt$Protein), " - Pept: ", unique(ddt$Peptide)))
  
  p.line <- ggplot(data=ddt, aes(x=Channel, y=Abundance, col=Run)) + 
    geom_point() + 
    geom_line(aes(group = Run )) +
    # scale_color_brewer(palette = "Dark2",na.value = "grey70") +
    facet_grid(PeptideN + Charge ~ Mixture, switch="y") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  ppt.p <- grid.arrange(p.tile, p.line, nrow = 2)
  return(ppt.p)
}
################################




################################
## plot proteins in tile and line
prt.plot <- function(ddt) {
  ddtl <- melt.data.table(ddt, id.vars = names(ddt[,-c("Abundance", "Abundance.noNorm")]), 
                          variable.name="Normalization", value.name="Abundance") #to long format
  ddtl[Normalization %in% c("Abundance", "Abundance.noNorm"), Normalization:= c("Norm.", "noNorm.")]    
  
  tmp <- data.table(Peptide = as.character(unique(ddt$Peptide)), PeptideN=NA)
  tmp[, PeptideN := sapply(Peptide, function(x) .spltpep(x,n=9))]
  ddt <- ddt[tmp, on="Peptide"]
  
  p.tile <- ggplot(data=ddt, aes(x=Run, y=PeptideN, fill=Abundance)) +
    geom_tile(color="white", size=0.05) +
    scale_fill_distiller(palette = "Spectral",na.value = "grey70") +
    facet_grid(PeptideN + Charge~Mixture+Channel, scales = "free", space = "free", switch="y") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.text.y=element_blank()) +
    ggtitle(paste0("Prot: ", unique(ddt$Protein) ))
  
  p.line <- ggplot(data=ddt, aes(x=Channel, y=Abundance, col=Run)) +
    geom_point() +
    geom_line(aes(group = Run )) +
    # scale_color_brewer(palette = "Dark2",na.value = "grey70") +
    facet_grid(Peptide + Charge ~ Mixture, switch="y") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  
  # ppt.p <- grid.arrange(p.tile, p.line, p.tile.n, p.line.n, nrow = 2)
  ppt.p <- grid.arrange(p.line, p.tile, nrow = 2)
  return(ppt.p)
} 
################################




################################
## plot protein profile 
plot.profile <- function(dd,ddd) { 
  tmp <- data.table(Peptide = as.character(unique(dd$Peptide)), PeptideN=NA)
  tmp[, PeptideN := sapply(Peptide, function(x) .spltpep(x,n=9))]
    dd <- dd[tmp, on="Peptide"]
    dd$PeptideSeq <- as.factor(dd$PeptideN)
      range <- round(summary(dd$Abundance))[c(1,6)]
  
      
  ggplot(dd, aes(x=Channel, y=Abundance, col=PeptideSeq)) + 
    
    geom_line(data=ddd, mapping=aes(x=Channel, y=Abundance, group=Gene), color="grey70", size=1.6) +
    geom_point(data=ddd, size=3, mapping=aes(x=Channel, y=Abundance), fill="grey70", color="black", stroke = 1) +
    
    geom_point(alpha=1) +
    geom_line(aes(group = Feature), alpha=1) +
    
    geom_text(aes(y=range[2]+1, label=Condition), col="grey20", hjust=1, vjust=-0.5, 
              check_overlap = TRUE,  cex=3, angle = 90) +
    geom_text(aes(y=range[2]+1, label=paste0("Rep.",BioReplicate)), col="grey20", 
              check_overlap = TRUE, hjust=1, vjust=1.5, cex=3, angle = 90) +
    
    scale_y_continuous(limits=c(range[1],range[2]+1 )) +
    
    facet_grid(. ~ factor(paste0("Mixture ",Mixture)), switch="y") +
    
    
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.position="right") + 
    
    ggtitle(paste0("Protein: ", unique(dd$Protein), "  -  Gene: ", unique(ddd$Gene)) ) +
    guides(fill=guide_legend(ncol=3, byrow=TRUE),
           colour = guide_legend(override.aes = list(alpha=1)))
}
################################




################################
## pca plot
pca.plot <- function(dc) {
  pcap <-
    fviz_pca_ind(dc,
                 col.ind = "cos2", # Color by the quality of representation
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE)
  return(pcap) }
################################




###############################
## plot heatmap
hmap.plot <- function(hm.mat0, reftb) {
  
  da <- data.table(cond=colnames(hm.mat0))
  da[, c("Condition", "BioReplicate") := tstrsplit(cond, " - Rep",  fixed=TRUE)]
  da[, c("Type", "Time") := tstrsplit(Condition, ".",  fixed=TRUE)]
  # da$Replicate <- gsub("Rep","", da$Replicate)
  da <- char2fact(da)
  reftb <- char2fact(reftb)
  da <- da[reftb[,c("Condition","BioReplicate", "Mixture")], on=c("Condition","BioReplicate")]
  
    daList <- as.list(NULL)
    daList[["Time + Type + BioReplicate"]] <- da[order(Time, Type,BioReplicate, decreasing=FALSE),]
    daList[["BioReplicate + Type + Time"]] <- da[order(BioReplicate,Type,Time, decreasing=FALSE),]
    daList[["Mixture + Type + Time"]] <- da[order(Mixture,Type,Time, decreasing=FALSE),]
    daList[["Type + Time"]] <- da[order(Type,Time, decreasing=FALSE),]
    
  
    hmList <- as.list(NULL)
    for (i in 1:length(daList) ) {  
      
      dda <- as.data.frame(daList[[i]])
                           
      hm.mat <- hm.mat0[,as.character(dda$cond)]
      nms <- data.table(names=rownames(hm.mat))
      nms[, col:="black"]
      nms[names %like% "Irf4" | names %like% "Nck1" | names %like% "Rela" | names %like% "Cblb" , col := "red"]
      dda <- as.data.frame(daList[[i]][,c("BioReplicate", "Mixture", "Type", "Time") ])
      
      ha = HeatmapAnnotation(df = dda, 
                             col = list(BioReplicate = c("3" = "thistle2", "7" = "tan1", "8" = "lightseagreen"),
                                        Mixture = c("2" = "honeydew2", "3" = "honeydew3", "4" = "steelblue2"),
                                        Time = c("0h" = "khaki", "24h" = "pink3", "48h" = "mediumpurple4"),
                                        Type = c("WT" = "steelblue4", "Cblb" = "tomato3")),
                             show_legend = c(TRUE, TRUE, TRUE, TRUE))
      
        
      if (i != length(daList)) {
        clustCol = FALSE
        rowName = FALSE 
        ttl = paste0("HM ordered by:\n", names(daList)[i])
        columnCol = rep(c("black","grey40"), nrow(dda)/2)
      } else {
        clustCol = TRUE
        rowName = TRUE   
        ttl = "HM ordered by\ncluster"
        columnCol = rep(c("black"), nrow(dda))
        }
      
      hmList[[names(daList)[i]]] <- Heatmap(hm.mat, name = "Abundance",
                         # column_dend_height = unit(2, "cm"), 
                         row_dend_width = unit(2, "cm"),
                         na_col = "white", 
                         clustering_distance_rows = "euclidean",
                         clustering_distance_columns = "euclidean",
                         cluster_columns = clustCol,
                         # col = circlize::colorRamp2(c(3,7), c("darkblue", "yellow"))
                         col <- colorRampPalette(rev(brewer.pal(10, "RdYlBu")) )(256),
                         row_names_gp = gpar(col = nms$col, fontsize = 6),
                         column_names_gp = gpar(col = columnCol, fontsize = 8),
                         column_title = ttl,
                         show_row_names = rowName,
                         column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                         top_annotation = ha) 
      
      
      }
  
  hm.all <- hmList[[1]] + hmList[[2]] + hmList[[3]] + hmList[[4]]
  
  return(hm.all)
}
##############################




##############################
## Volcano plots: all in one page
volc.p.all <- function(ldt) {
  
    colvec <- c("#899DA4", "#efe6c2", "#ffa75b" ,"#C93312")
    ldt[, col := colvec[1]]
    ldt[ adj.P.Val < 0.05 , col := colvec[2]]
    ldt[ abs(logFC) > 1   , col := colvec[3]]
    ldt[ adj.P.Val < 0.05 & abs(logFC) > 1 , col := colvec[4]]
    ldt$col <- as.factor(ldt$col)
    
    tmp <- ldt[ adj.P.Val < 0.05, .SD[which.max(adj.P.Val)], by = Comparison]
    tmp[, cutoff := -log10(P.Value)]
    ldt <- merge(ldt, tmp[, c("Comparison", "cutoff")], all.x = TRUE)
    
    range <- c( floor(summary(fitList$logFC)["Min."]), ceiling(summary(fitList$logFC)["Max."]) )
    ticks <- c(-8,-4,-2,-1,0,1,2,4,8)
    
    cntb <- ldt[, .N, by=list(Comparison, col)]
    cntb[, perc := round(100*N/sum(N),2 ), by=Comparison]
    
        maxY <- round(max(-log10(ldt$P.Value), na.rm=TRUE) )
        maxX <-   as.numeric(range[-1])
        cntb[, xx := maxX]
        cntb[, yy := seq(col), by=Comparison ]
        cntb$yy <- 0.7*as.numeric(cntb$yy)
        
    # FC = paste("LogFC > |1|", ), 
    # P.adj=paste("P.adj<", AdjustedCutoff, sep=""), 
    # FC_P.adj=paste("P.adj<", AdjustedCutoff, " & LogFC>|", FCCutoff, "|", sep=""))) +
          
    ## order facets
    # ldt$facet <- factor(as.factor(ldt$Comparison), levels = c(
    #                                     "WT-Cblb"          ,"Cblb.0h-WT.0h"  , "Cblb.24h-WT.0h" , "Cblb.48h-WT.0h" , 
    #                                     "Cblb.0h-Cblb.24h" , "Cblb.0h-WT.24h", "Cblb.24h-WT.24h", "Cblb.48h-WT.24h",
    #                                     "Cblb.0h-Cblb.48h" , "Cblb.0h-WT.48h", "Cblb.24h-WT.48h", "Cblb.48h-WT.48h",
    #                                     "Cblb.24h-Cblb.48h", "WT.0h-WT.24h"  , "WT.0h-WT.48h"   , "WT.24h-WT.48h" ))
    # 
    # ldt  <- ldt[order(match(Comparison, as.factor(facet)))]
    # cntb <- cntb[order(match(Comparison, as.factor(facet)))]
    # setnames(cntb, new = "facet", old="Comparison")
    
    vplot <- 
      ggplot(data=ldt, aes(x=logFC, y = -log10(P.Value), colour = col, group = col)) +
      geom_point(alpha=0.6, size=1.75) +
      xlab("log2 fold change") + ylab("-log10 p-value") +
      # scale_color_manual(values=wes_palette(n=4, name="Royal1"))
      # scale_color_manual(values=c("#899DA4", "#efe6c2", "#ffa75b" ,"#C93312")) +
      scale_color_manual(values=levels(ldt$col)) +
      facet_wrap(~Comparison, ncol=4) +
    
      guides(colour = guide_legend(override.aes=list(size=2.5))) +
      
      xlab(bquote(~Log[2]~ "fold change")) +
      ylab(bquote(~-Log[10]~italic(P))) +
      
      geom_vline(xintercept=c(-1, 1), linetype="dashed", colour="grey40", size=0.4) +
      geom_vline(xintercept=c(-2, 2), linetype="longdash", colour="grey70", size=0.4) +
      geom_vline(xintercept=0, linetype="dotted", colour="grey", size=0.4) +
      geom_hline(aes(yintercept=cutoff, group = Comparison), linetype="dashed", colour="grey40", size=0.4) +
      
      scale_x_continuous(breaks = ticks, limits=range) +
    
      theme_bw(base_size=24) +
      theme(legend.background=element_rect(),
            plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
            
            panel.grid.major=element_blank(),	#Remove gridlines
            panel.grid.minor=element_blank(),	#Remove gridlines
            
            axis.text.x=element_text(angle=0, size=10, vjust=1),
            axis.text.y=element_text(angle=0, size=10, vjust=1),
            axis.title=element_text(size=10),
            strip.text.x = element_text(size=10, margin = margin(1,0,1,0)),
            
            legend.position="none",		
            legend.key=element_blank(),		
            legend.text=element_text(size=8),	
            title=element_text(size=8))

  sub <- ldt[Gene %in% c("Irf4" , "Nck1" , "Rela" , "Cblb" ), ]
  sub$side <- ifelse(sub$logFC < 0, -1, 1)
  
  vplot <- vplot + 
  geom_point(data= sub, size=3, mapping=aes(x=logFC, y = -log10(P.Value)), 
                   color="black", stroke = 1, shape=0) +
    
  geom_label_repel(data=sub, aes(label = Gene),
                   force = 7,
                   show.legend = FALSE,
                   label.size = 0.1,
                   alpha = 0.7,
                   # nudge_x = ifelse(sub$side == 1, 4, -4),
                   point.padding = unit(0.5, "lines"),
                   box.padding = unit(0.1, "lines"),
                   segment.alpha =1 ,
                   segment.colour = "grey30") +
  geom_label(data=cntb, aes(x=xx, y=yy, label = paste("#", N), col= col), hjust=1, size=3) 
  
  return(vplot)
}




##############################
## Volcano plots: in seperate pages
volc.p.ind <- function(fitList, topList) {
  
  ldt <- copy(fitList)
  top <- copy(topList)

  colvec <- c("#899DA4", "#efe6c2", "#ffa75b" ,"#C93312")
  ldt[, col := colvec[1]]
  ldt[ adj.P.Val < 0.05 , col := colvec[2]]
  ldt[ abs(logFC) > 1   , col := colvec[3]]
  ldt[ adj.P.Val < 0.05 & abs(logFC) > 1 , col := colvec[4]]
  ldt$col <- as.factor(ldt$col)
  
  tmp <- ldt[ adj.P.Val < 0.05, .SD[which.max(adj.P.Val)], by = Comparison]
  # tmp <- ldt[ adj.P.Val < 0.05, .SD[which.max(-log10(P.Value))], by = Comparison]
  tmp[, cutoff := -log10(P.Value)]
  ldt <- merge(ldt, tmp[, c("Comparison", "cutoff")], all.x = TRUE)
  
  range <- c( floor(summary(fitList$logFC)["Min."]), ceiling(summary(fitList$logFC)["Max."]) )
  ticks <- c(-8,-4,-2,-1,0,1,2,4,8)
  
  
  cntb <- ldt[, .N, by=list(Comparison, col)]
  cntb[, perc := round(100*N/sum(N),2 ), by=Comparison]

    
  volcList <- list()
  
  for (ix in 1:uniqueN(ldt$Comparison)) {
    
    subldt <- ldt[Comparison == unique(ldt$Comparison)[ix] ]
    subtop <- top[Comparison == unique(ldt$Comparison)[ix] & id %in% 1:30]
    
    subcntb <- cntb[Comparison == unique(cntb$Comparison)[ix] ]
      maxY <- round(max(-log10(subldt$P.Value), na.rm=TRUE) )
      maxX <-   as.numeric(range[-1])
        subcntb[, xx := maxX]
        subcntb[, yy := 0.4*(seq(nrow(subcntb))) ]
    
      volcList[[ix]]  <-
        ggplot(data=subldt, aes(x=logFC, y = -log10(P.Value), colour = col, label=Gene)) +
        geom_point(alpha=0.5, size=1.75) +
        xlab("log2 fold change") + ylab("-log10 p-value") +
        # scale_color_manual(values=wes_palette(n=4, name="Royal1"))
        # scale_color_manual(values=c("#899DA4", "#efe6c2", "#ffa75b" ,"#C93312")) +
        # scale_color_manual(values=unique(subldt$col)) +
        scale_colour_identity() +
        
        geom_label(data=subcntb, aes(x=xx, y=yy, label = paste("#", N), col= col), hjust=1, size=5) +
        
        
        guides(colour = guide_legend(override.aes=list(size=2.5))) +
        
        xlab(bquote(~Log[2]~ "fold change")) +
        ylab(bquote(~-Log[10]~italic(P))) +
        
        geom_vline(xintercept=c(-1, 1), linetype="dashed", colour="grey40", size=0.4) +
        geom_vline(xintercept=c(-2, 2), linetype="longdash", colour="grey70", size=0.4) +
        geom_vline(xintercept=0, linetype="dotted", colour="grey", size=0.4) +
        geom_hline(aes(yintercept=cutoff, group = Comparison), linetype="dashed", colour="grey40", size=0.4) +
        
        scale_x_continuous(breaks = ticks, limits=range) +
        
        theme_bw(base_size=24) +
        theme(legend.background=element_rect(),
              plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
              
              panel.grid.major=element_blank(),	#Remove gridlines
              panel.grid.minor=element_blank(),	#Remove gridlines
              
              axis.text.x=element_text(angle=0, size=10, vjust=1),
              axis.text.y=element_text(angle=0, size=10, vjust=1),
              axis.title=element_text(size=10),
              strip.text.x = element_text(size=10, margin = margin(1,0,1,0)),
              
              legend.position="none",		
              legend.key=element_blank(),		
              legend.text=element_text(size=8),	
              title=element_text(size=8)) +
        
        ggtitle(unique(subldt$Comparison)) 
      
      
      # sub <- subldt[adj.P.Val < 0.05 & abs(logFC) > 1 , ]
      # sub$side <- ifelse(sub$logFC < 0, -1, 1)
      subsub <- subldt[Protein %in% subtop$Protein, ]
      subsub$side <- ifelse(subsub$logFC < 0, -1, 1)
      
        
        volcList[[ix]] <-  volcList[[ix]] +
        geom_point(data= subsub, size=3, mapping=aes(x=logFC, y = -log10(P.Value)), 
                   color="black", stroke = 1, shape=0) +
        # facet_wrap(~Comparison, ncol=4, as.table=TRUE) +
        geom_label_repel(data=subsub, aes(label = Gene),
                         force = 7,
                         show.legend = FALSE,
                         label.size = 0.2,
                         alpha = 1,
                         # nudge_x = ifelse(sub$side == 1, 4, -4),
                         point.padding = unit(0.5, "lines"),
                         box.padding = unit(0.1, "lines"),
                         segment.alpha =1 ,
                         segment.colour = "grey30") 
        
        
  }
  return(volcList)
}




pBox <- function(wt, x, y, col, ...) {
  gp <- 
    ggplot(wt, aes_string(x=x, y=y, col=col)) + 
    geom_boxplot() +
    scale_color_brewer(palette = "Dark2") + 
    theme_bw() + 
    facet_wrap(~Mixture, scales = "free_x") +
    # theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position="none") 
  
  ggtitle(y)
  return(gp)
}


pDen <- function(wt, x, col, ...) {
  gp <- 
    ggplot(wt, aes_string(x=x, col=col)) + 
    geom_density(aes(linetype = Channel)) +    
    # scale_color_brewer() + 
    theme_bw() + 
    facet_wrap(~Mixture) +
    theme(legend.position = "none") +
    ggtitle(x) + 
    xlim(0,12)
  return(gp)
}












