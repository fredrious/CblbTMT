

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
hmap.plot <- function(hm.mat) {
  require("RColorBrewer")
  
  da <- data.table(cond=colnames(hm.mat))
  da[, c("Condition", "Replicate") := tstrsplit(cond, "\n",  fixed=TRUE)]
  da[, c("Type", "Time") := tstrsplit(Condition, ".",  fixed=TRUE)]
  da$Replicate <- gsub("Rep","", da$Replicate)
  da <- da[reftb[,c("Condition","Replicate", "Mixture")], on=c("Condition","Replicate")]
  
  ha = HeatmapAnnotation(df = as.data.frame(da[,c("Replicate", "Mixture", "Type") ]), 
                         col = list(Replicate = c("3" = "blue", "7" = "grey", "8" = "red"),
                                    Mixture = c("2" = "green", "3" = "yellow", "4" = "pink"),
                                    Type = c("WT" = "darkblue", "Cblb" = "darkred")),
                         show_legend = c(TRUE, TRUE, TRUE))
  
  hm.plot <- Heatmap(hm.mat, name = "Abundance", 
                     # column_dend_height = unit(2, "cm"), 
                     row_dend_width = unit(2, "cm"),
                     na_col = "white", 
                     clustering_distance_rows = "pearson",
                     cluster_columns = FALSE,
                     # col = circlize::colorRamp2(c(3,7), c("darkblue", "yellow"))
                     col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256),
                     top_annotation = ha)
  
  hmallclust.plot <- Heatmap(hm.mat, name = "Abundance", 
                             # column_dend_height = unit(2, "cm"), 
                             row_dend_width = unit(2, "cm"),
                             na_col = "white", 
                             clustering_distance_rows = "pearson",
                             cluster_columns = TRUE,
                             # col = circlize::colorRamp2(c(3,7), c("darkblue", "yellow"))
                             col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256),
                             top_annotation = ha)
  
  return(list(hm.plot, hmallclust.plot))
}
##############################


