

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
  ggplot(dd, aes(x=Channel, y=Abundance, col=factor(PeptideN))) + 
    
    geom_line(data=ddd, mapping=aes(x=Channel, y=Abundance, group=MPmethod), 
              color="grey70", size=1.6) +
    geom_point(data=ddd, size=3, mapping=aes(x=Channel, y=Abundance, shape=MPmethod), 
               fill="grey70", color="black", stroke = 1) +
    
    scale_shape_manual(name = "MP method", values = c(22,24),
                       breaks = c("allMix", "indMix"),
                       labels = c("all Mix", "each Mix")) +
    
    geom_point(alpha=1) +
    geom_line(aes(group = Feature), alpha=1) +
    
    geom_text(aes(y=max(Abundance), label=Condition), hjust=1, vjust=-0.5, 
              check_overlap = TRUE, col="black", cex=3, angle = 90) +
    geom_text(aes(y=max(Abundance), label=paste0("Rep.",BioReplicate)), 
              check_overlap = TRUE, hjust=1, vjust=1.5, col="black", cex=3, angle = 90) +
    
    # scale_colour_brewer(type = "qual", palette = "Dark2") +
    
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
    fviz_pca_ind(pca2,
                 col.ind = "cos2", # Color by the quality of representation
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE)
  return(pcap) }
################################




###############################
## plot heatmap
hmap.plot <- function(hm.mat0) {
  require("RColorBrewer")
  
  da <- data.table(cond=colnames(hm.mat0))
  da[, c("Condition", "BioReplicate") := tstrsplit(cond, "\nRep",  fixed=TRUE)]
  da[, c("Type", "Time") := tstrsplit(Condition, ".",  fixed=TRUE)]
  # da$Replicate <- gsub("Rep","", da$Replicate)
  da <- da[reftb[,c("Condition","BioReplicate", "Mixture")], on=c("Condition","BioReplicate")]
  da <- char2fact(da)
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
        ttl = paste0("HM ordered by:   ", names(daList)[i])
        columnCol = rep(c("black","grey40"), nrow(dda)/2)
      } else {
        clustCol = TRUE
        rowName = TRUE   
        ttl = "HM ordered by cluster"
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
                         row_names_gp = gpar(col = nms$col, fontsize = 8),
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


