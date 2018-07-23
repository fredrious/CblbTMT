


rm(list=ls())

ddir <- "/Users/farhad/_Rspace/cblb_wtVkn/data/"
rdir <- "/Users/farhad/_Rspace/cblb_wtVkn/results/"
rdir <- "/Users/farhad/_Rspace/cblb_wtVkn/scripts/"
setwd("/Users/farhad/_Rspace/cblb_wtVkn/newMeasurements/")

library("data.table")
library("openxlsx")
library("dplyr")
library("tidyr")
library("ggplot2")
library("limma")
library("lme4")
library("vsn")
library("wesanderson")
library("readxl")
library("gridExtra")
library("biomaRt")
library("hpar")
library("ComplexHeatmap")
library("circlize")

library("devtools")
load_all(paste0("/Users/farhad/_Rspace/_loclib/MSstatsTMT_master/R"))



################################### slef-defined functions
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


## Function to split long peptides sequences into smaller pieces of length n (by default n=8)
## and insert \n at split location. This is for having shorter texts on plots. In practice, 
## \n breaks peptide names into multiple lines.
.spltpep <- function(x, n=8) {
  nc <- nchar(x)
    if (nc > n) {
    pos <- seq(nc)[(seq(nc) %% n) == 0]
      if (length(pos) > 1) {
        if (pos[-1] == nc) pos <- head(pos,-1)
      }
    from <- c(1,pos+1)
    to <- c(pos, nchar(x))
    xx <- substring(x, from, to)
    st <- do.call(paste, c(as.list(xx), sep = "\n")) 
    } else { st = x }
    return(st)
  }



## function: median polish
f.medpol <-  function(dw, para) {
  clnm <- colnames(dw)
  dw <- as.matrix(dw)
  mp  <-  stats::medpolish(dw, na.rm=TRUE, trace.iter = FALSE)
  if (toupper(para)=="MEDRES") {
    tmp <- mp$overall + apply(mp$residuals, 2, function(x) median(x, na.rm = TRUE)) # Abundance = overall median + median(residual)
  } else { tmp <- mp$overall + mp$col } # Abundance = overall median + column effect
  result <- data.table(Channel=clnm, Abundance=tmp)
  res <- mp$residual
  return(list(result,res))
}



## Function to get the median value
which.medX <-  function(x) {
        xx <-  x >= median(x, na.rm = TRUE)
        which( x == x[xx][which.min(x[xx])] )}
## gives the median if length(x)=2*n+1 or the larger of the two middle values if length(x)=2*n. 


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

 ## ensemble mouse
 ens <- useMart("ensembl", "mmusculus_gene_ensembl")

 
 
 
 ###################################
 ## Read in PSM data
 meas <- "new"
 if (meas=="new") {
   ## Reading PSM data --> new measurements. 25.06.2018 --> run 39 is missing!
   psm.raw <- as.data.table(read_excel("180504_P_268_IC_Pool3_7_8_SPS_180620_PSMs.xlsx"))
   ## Rename Spectrum files to match format
   psm.raw$`Spectrum File` <- gsub("_SPS","", psm.raw$`Spectrum File`)
   psm.raw$`Spectrum File` <- gsub("_tr4","", psm.raw$`Spectrum File`)
 } else {
   ## Reading PSM data
   psm.raw <- as.data.table(
     read.xlsx(paste0(ddir,"cblb_raw/171115_P_268_IC_allPools_PSMs.xlsx"),
               sheet=1, startRow=1, colNames=TRUE,
               rowNames=FALSE, skipEmptyRows=FALSE)) 
   ## Removing replicate 2 --> File.ID in [F1.1 - F1.12]
   psm.raw <- psm.raw[!psm.raw$File.ID %like% "F1.", ] 
   
 }
 
 psm.tmp <- copy(psm.raw)
 names(psm.tmp) <- gsub(" ",".", names(psm.tmp))
 psm.tmp[, Run:= lapply(.SD, function(x) as.integer(sub(".*\\_([^.]+)\\..*", "\\1", x))), .SDcols="Spectrum.File"]
 i.x <- which(names(psm.tmp) %like% "Abundance")
 setnames(psm.tmp, old = names(psm.tmp)[i.x], new = gsub("Abundance:.", "", names(psm.tmp)[i.x]))
 psm.dt <- melt.data.table(psm.tmp, id.vars = names(psm.tmp)[!names(psm.tmp) %in% c("126", "127", "128", "129", "130", "131")],
                           variable.name="Channel", 
                           value.name="Intensity")
 ###################################
 
 
 
 ###################################
 ## Reading PSM Annotation table and reformatting 
 psm.ano <- fread(file=paste0(ddir,"cblb_raw/171115_P_268_IC_Annotation.csv"), sep=c(";"), header=TRUE, dec=".", key = "RAW_Filenummer")
 psm.ano[, c("from", "to") := tstrsplit(RAW_Filenummer, "-",  fixed=TRUE)]
 psm.ano[, c("from","to") := lapply(.SD, as.numeric), .SDcols=c("from","to")]
 # psm.ano[, Run.i:= lapply(.SD, function(x) as.numeric(sub(".*\\_([^.]+)\\..*", "\\1", x))), .SDcols=Run]
 psm.ano <- psm.ano[, list(Run = from:to), by = .(RAW_Filenummer, Channel, Condition, Mixture, BioReplicate)]
 psm.ano[, RAW_Filenummer := NULL]
 ###################################
 
 
 
 ###################################
 ## change Charachter and Integer columns to Factor
 ## PSM table:
 changeCols<- c(names(Filter(is.character, psm.dt)), names(Filter(is.integer, psm.dt)))
 psm.dt[,(changeCols):=lapply(.SD, as.factor),.SDcols=changeCols]
 ## Annotation table:
 changeCols<- c(names(Filter(is.character, psm.ano)), names(Filter(is.integer, psm.ano)))
 psm.ano[,(changeCols):=lapply(.SD, as.factor),.SDcols=changeCols]
 ###################################
 
 
 
 
## Dynamic parameters for version chack
###################################
# meas.arr <- c("old" ,"new")                   ## Which measurement
PSM2PPT.arr = c(TRUE, FALSE)                  ## PSM to Peptide (Ratio > Median) 
FractComb.arr <-  c("med", "max", "none")      ## combine fractions
remMix2.arr <- c(TRUE, FALSE)                 ## remove mixture 2 or not (high missingness in mix2)
allChannelsIn.arr <- c(TRUE, FALSE)           ## whether or not consider Features with complete Channel set
MPpara.arr <-  c("coleff","medres")           ## medianpolish:  MPpara=coleff --> Abundance = overall median + column effect
                                              ## medianpolish:  MPpara=medres --> Abundance = overall median + median(residual)
medpolON.arr <-  c("allMix", "indMix")        ## medianpolish to be applied to whole protein block (allMix), or to each Mixture seperately (indMix)
whichPSM.arr <-c("med", "max")                ## for the same PSM per run, mixture and channel
TopNperc <-  5                                ## Limma TopTable: top N percent hits

onlySharedPRT.arr <- c(TRUE, FALSE)           ## Reduce dataset to only shared Proteins between all Mixtures
minPptMIX.arr <- c(1,2,0)                       ## minimum number of Peptide in each Mixture
mix3rev.arr <- c(TRUE, FALSE)                    




## Constatnt parameters for version chack
PSM2PPT <- FALSE
MPpara <-  "coleff"
medpolON <-  "allMix"
allChannelsIn <- FALSE
whichPSM <- "med"
remMix2 <- TRUE      
minPptMIX <- 1
FractComb <- "max"
onlySharedPRT <- FALSE

mix3rev <- TRUE


# for (ix in 1:length(meas.arr)) {
  # for (jx in 1:length(FractComb.arr)) {
    # for (kx in 1:length(remMix2.arr)) {
      for (lx in 1:length(allChannelsIn.arr)) {
        # for (mx in 1:length(medpolON.arr)) {
          # for (nx in 1:length(whichPSM.arr)) {
            # for (px in 1:length(onlySharedPRT.arr)) {
              for (qx in 1:length(minPptMIX.arr)) {
                
                # meas <- meas[ix]
                # FractComb <- FractComb.arr[jx]
                # remMix2 <- remMix2.arr[kx]
                allChannelsIn <- allChannelsIn.arr[lx]
                # medpolON <- medpolON[mx]
                # whichPSM <- whichPSM.arr[nx]
                # onlySharedPRT <- onlySharedPRT.arr[px]
                minPptMIX <- minPptMIX.arr[qx]
                

                
                
########################################################
################################### Filterring procedure
  ## Merge with annotation
  work0 <- merge(psm.dt, psm.ano, by = c("Run","Channel"), allow.cartesian=TRUE)
  work0 <-  work0[, c("Run", "Channel", "Master.Protein.Accessions", "#.Protein.Groups", "Annotated.Sequence", "Charge", 
                      "Ions.Score", "Rank", "MS.Order", "File.ID", "RT.[min]", "Condition", "Quan.Info",
                      "Mixture", "BioReplicate", "Intensity")]
  setnames(work0, old = c("Master.Protein.Accessions", "Annotated.Sequence", "#.Protein.Groups"), 
           new = c("Protein", "Peptide", "numProtein") )
  work0[Intensity < 1, Intensity := 1]
                
  ####### STEP 1: Removing shared proteins
  ## Number of protein accessions in Protein --> this should be equal to 1.
  work0[, cntProt := lengths(strsplit(as.character(Protein), ";"))]
  # Remove rows with Protein == "sp"
  work0 <- work0[!Protein=="sp", ] 
  # Remove rows with more than 1 protein accession and no "sp" (shared protein)
  work0 <- work0[cntProt>1 & Protein %like% "sp" | cntProt==1, ]
  ## For Protein containing "sp" and one protein accession, remove "sp" and stay with protein accession.
  ## If the number of Protein accessions after removing "sp" is more than 1 --> shared protein --> remove row. 
  work0[,Protein := gsub("sp;","", Protein, fixed = TRUE)]
  work0[,Protein := gsub("; sp","", Protein, fixed = TRUE)]
  work0[,Protein := gsub(" ","", Protein, fixed = TRUE)]
  work0 <- work0[!Protein=="sp", ] 
  
  work0[, cntProt := lengths(strsplit(as.character(Protein), ";"))]
  work0 <- work0[cntProt==1]
  ## The process remove 150462 rows --> 11% of the data
  # work0[,Quan.Info:=factor("Unique")] #Now Quant.info can be marked as unique. ## this id needed for MSstatsTMT
  
  work0[,cntProt:=NULL]
  ## define Feature and Proptide
  work0[, Feature:= do.call(paste, c(.SD, sep = "_")), .SDcols=c("Protein", "Peptide", "Charge")]
  work0[, Proptide:= do.call(paste, c(.SD, sep = "_")), .SDcols=c("Protein", "Peptide")]
  ## Define Abundance as log2(Intensity)
  work0[, Abundance := log2(Intensity)]
  # Continue with Abundance: Andreas 6.6.2018
  work0 <- work0[, list(Run, Channel, Protein, Peptide, Feature, 
                       Charge, Proptide, Condition, BioReplicate, Mixture, Abundance)]
  
  
  ####### STEP 2: removing PSMs with protein-counts != 1
  # work <- work[numProtein == 1, ] # This is not neccessary after removing shared proteins manually in STEP 1.
  # work0[,numProtein := NULL] # After STEP 1, protein counts per PSM are already 1.
  
  
  ####### STEP 3: Peptides, that are used in more than one proteins (Similar to previous step. This is a double check!)
  notUnqPpt <- work0[, c(.(cnt = uniqueN(Protein))), by = list(Peptide)] [cnt!=1, Peptide]
  if (length(notUnqPpt) > 0) work0 <- work0[!Peptide %chin% notUnqPpt,]


  
  work1 <- copy(work0)
  work1[Abundance=="NaN", Abundance := NA]
  
  
  ## NA count per feature and Run
  NAcnt.work1 <- work1[, c(.(NAper = sapply( sum(is.na(Abundance)) / length(Abundance) , function(x) round(x,3) ) ) ,
                           .(NAcnt = sum(is.na(Abundance))),
                           .(TTLcnt = length(Abundance))) ,
                         by = list(Feature )]
  # table(NAcnt.work1$NAcnt)
  # table(NAcnt.work1$TTLcnt)
  NA.tb <- as.data.table(table(NAcnt.work1$NAper))
  NAfeat <- round(100 * (sum(NA.tb[V1!=0, N]) / sum(NA.tb$N)) ,3) # percentage of Features witj at least 1 NA
  mnNA <- round(100 * (summary(as.numeric(NA.tb[V1!=0, V1]))["Min."]),4)
  mxNA <- round(100 * (summary(as.numeric(NA.tb[V1!=0, V1]))["Max."]),4)
  
  print(paste0(NAfeat, "% of the Features have at least 1 NA value across all Runs and Channels."))
  print(paste0("Degree of missingness in Features varies from ~",mnNA, "% to ~", mxNA, "%"))
  
  
  ########################################################
  ################################### Dealing with multiple PSM measurements within and across fractions.

  ## Below: Number of measurements (intensity values) per Feature per Run. This must ideally be equal to number of Channels (here 6)
  ## We have to decide what to do when a Feature is measured more than once. 
  cnt.tmp <-  work1[, c(.(cnt = length(Abundance))), by = c("Run","Feature")] 
  ## table(cnt.tmp$cnt) --> contingency table showing row-counts at each count level 

      
      ###############################
      ####### STEP 1: mult. measurements for same Feature, same Mixture, same Channel and same Run: calc. Mean OR Median
      ###############################
      if (remMix2) { work1 <-  work1[Mixture != 2] } # remove mixture 2
  # work11 <- copy(work1)
      
      if (toupper(whichPSM) == "MEAN") {
        work1[, Abundance := mean(Abundance, na.rm=TRUE), by=list(Feature, Mixture, Run, Channel)] #Mean
      }
      if (toupper(whichPSM) == "MED") {
        work1[ ,Abundance := median(Abundance, na.rm=TRUE), by=list(Feature, Mixture, Run, Channel)] #Median
      }
      if (toupper(whichPSM) == "MAX") {        
        work1[ ,Abundance := max(Abundance, na.rm=TRUE), by=list(Feature, Mixture, Run, Channel)] #max
        # work1[work1[, .I[which.max(Abundance)],  by=list(Feature, Mixture, Run, Channel)]$V1]
      } 
        ## after mean/max/median ... work is not unique:
        # work1[Feature == "P17742_sIYGEKFEDENFILk_3" & Run == 40 & Channel ==126,]
        # Run Channel Protein         Peptide                  Feature Charge               Proptide
        # 1:  40     126  P17742 sIYGEKFEDENFILk P17742_sIYGEKFEDENFILk_3      3 P17742_sIYGEKFEDENFILk
        # 2:  40     126  P17742 sIYGEKFEDENFILk P17742_sIYGEKFEDENFILk_3      3 P17742_sIYGEKFEDENFILk
        # Condition BioReplicate Mixture Abundance
        # 1:     WT.0h            3       2  9.176075
        # 2:     WT.0h            3       2  9.176075
      
      work1 <- unique(work1)
      ## Example 1. BEFORE AVERAGING
      ## work0[Run==100 & Mixture==4 & Channel==130 & Feature =="P62301_vLPPNWk_2"]
      ## Example 1. AFTER AVERAGING
      ## work[Run==100 & Mixture==4 & Channel==130 & Feature =="P62301_vLPPNWk_2"]
      
      ## Number of Charge states per Peptide per Protein
      cnt.work1 <- work1[, .(n = uniqueN(Charge)), by = list(Proptide)] 
      # table(cnt.work1$n)
      #     1     2     3     4
      # 49961  4439   107     1
            ## @Andrea: Below examples are interesting in terms of imputation!
            ## I would suggest to remove PSM with a given percentage of missingness, like 25% or 50% !!
            ## Otherwise imputation may apply more error, than useful information. We can discuss.
            # prtp4 <- cnt.work1[n == 4, Proptide][1] # peptide with 4 Charge states
            # dt4 <- work1[Proptide == prtp4, ] # corresponding data block
            # 
            # prtp3 <- cnt.work1[n == 3, Proptide][1] # peptide with 3 Charge states
            # dt3 <- work0[Proptide == prtp3, ] # corresponding data block
            # 
            # prtp2 <- cnt.work1[n == 2, Proptide][1] # peptide with 2 Charge states
            # dt2 <- work0[Proptide == prtp2, ] # corresponding data block
            # 
            # ppt.plot(dt2)
            # ppt.plot(dt3)
            # ppt.plot(dt4)
            # 
            # dt <- work0[Protein == "Q3TTA7"]
            # ppt.plot(dt)
            
            
            
            # ##  # peptide with 2 Charge states
            # prtp2 <- cnt.work1[n == 2, Proptide][1] # peptide with 2 Charge states
            # dt2 <- work0[Proptide == prtp2 & Run == 91, ]
            # 
            # dt2[, Abundance1 := mean(Abundance, na.rm=TRUE), by=list(Feature, Mixture, Run, Channel)] #Mean
            # dt2[ ,Abundance2 := median(Abundance, na.rm=TRUE), by=list(Feature, Mixture, Run, Channel)] #Median
            # dt2[ ,Abundance3 := max(Abundance, na.rm=TRUE), by=list(Feature, Mixture, Run, Channel)] #max
            # 
            # 
            # ggplot(data=dt2, aes(x=Channel, col=Run)) + 
            #   geom_point(aes(y=Abundance)) + 
            #   geom_point(aes(y=Abundance1), col="darkred", size=3, shape=0) +
            #   geom_point(aes(y=Abundance2), col="black", size=3, shape=0) +
            #   geom_point(aes(y=Abundance3), col="blue", size=3, shape=0) +
            #   
            #   geom_line(aes(y=Abundance, group = Feature )) +
            #   # scale_color_brewer(palette = "Dark2",na.value = "grey70") +
            #   facet_grid(Peptide + Charge ~ Mixture, switch="y") +
            #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
            
            
            
            
            
      
      ###############################
      ####### STEP 2: From PSM to Peptide level:
      ###############################
      if (PSM2PPT) {
          ## CASE 1: Peptides with 2 PSM > i. compute PSM ratios (in log scale: ratio = PSM1 - PSM2)
          ##                               ii. compute median of ratios (medRatio)
          ##                               iii. when PSM2 is NA & PSM1 nonNA: subtract medRatio from PSM1 (PSM2 = PSM1 - medRatio)
          ##                               iv . when PSM1 is NA & PSM2 nonNA: add medRatio to PSM2 (PSM1 = PSM2 + medRatio)
          
        
          ## function to i. calculate PSM ratios
          ##             ii. impute NA for PSM within Peptide, Mixture, Run and Channel
          ##             iii. return original Abundance, imputed Abundance & median Abindance
          ## note: median abundance inserted in datatable for both charge states, meaning median abundance is duplicated for each peptide.
          ## This is useful for ggplot to have all three levels of Abundance in the same datatable.
          f.psm2ppt <-  function(ww) {
            ww <- t(ww[,"Feature"])
            colnames(ww) <- t(ww[,"Feature"])
            ww[is.na(ww[,2]),2] <- ww[is.na(ww[,2]),1] - median(ww[,1]-ww[,2], na.rm=TRUE) #filling NAs of 2nd Charge by nonNAs of 1st Charge 
            ww[is.na(ww[,1]),1] <- ww[is.na(ww[,1]),2] + median(ww[,1]-ww[,2], na.rm=TRUE) #filling NAs of 1st Charge by nonNAs of 2nd Charge 
            ww.med <- data.table(nm=rownames(ww), ww, medPeptide=apply(ww, 1, function(x) median(x, na.rm=TRUE))) #table of imputed and mean PSM
        
            PSM.imp <- melt.data.table( ww.med[,-4], id.vars = "nm", variable.name="Feature", value.name="impAbundance") #to wide format
            PSM.imp[, c("Mixture", "Run", "Channel") := tstrsplit(nm, "_", fixed=TRUE)] # table of imputed abundance
            PSM.imp[, c("Protein", "Peptide", "Charge") := tstrsplit(Feature, "_", fixed=TRUE)]
            
            PSM.med <- melt.data.table( ww.med[,c(1,4)], id.vars = "nm", value.name="medAbundance")
            PPT.dt <- PSM.imp[PSM.med[,-"variable"], on=c("nm")] # add median abundance
            
            return(PPT.dt)
          }
         
          ## subset data to Peptides with 2 PSMs
          pptpsm2 <- cnt.work1[n == 2, Proptide] # Peptides with 2 PSMs
          work.2psm <- work1[Proptide %in% pptpsm2, ] # corresponding data block
          work.2psm.w <- data.table::dcast(work.2psm, Proptide+Feature ~ Mixture + Run + Channel, value.var = "Abundance") #to wide format
          work.2psm.merged <- work.2psm.w[, f.psm2ppt(.SD) , by = .(Proptide)  ]  
          
          ## merge with original data and output
          work1 <- merge(work.2psm.merged, work.2psm[,c("Feature", "Mixture", "Run", "Channel", "Abundance")],
                       by = c("Feature", "Mixture", "Run", "Channel"), all.x = TRUE)
          
      } # end of PSM2PPT
     
      
      
  ########################################################
  ################################### Fraction combination: Median, Sum, Max, NONE!(whichPSM == "median")
  if (toupper(FractComb) != "NONE") {
    ##!!! UPDATE 14.06.2015: Marc suggested that instead of removing overlapped PSM between Runs of
    ##!!! the same biological mixture, take Maximum of the redundant values.
    work2 <- copy(work1)
    sharedPSM.Mix <- work2[, c(.(cnt = uniqueN(Run))), by = list(Feature, Mixture)] #table(sharedPSM.Mix$cnt)
    # check: ppt.plot(work[Proptide == "P52480_lAPITSDPTEAAAVGAVEASFk"])
    wch0 <- work2[sharedPSM.Mix[cnt==1,], on=c("Feature", "Mixture")] #unique measurements
    wch  <- work2[sharedPSM.Mix[cnt>1,], on=c("Feature", "Mixture")]  #redundant measurements
      wch1 <- wch[ , .(mean_int = mean(Abundance, na.rm=TRUE)), by = c("Protein","Feature","Mixture","Run")] # mean Channels for each PSM, Run and Mix
      
      if (toupper(FractComb) == "MAX") { #get the maximum fraction
          wch.frc  <- wch1[wch1[, .I[which.max(mean_int)],  by=c("Protein","Feature","Mixture")]$V1]}
      if (toupper(FractComb) == "MED") { #get the median fraction > for 2n fractions, get the 1st fraction above median > see: which.medX
          wch.frc <- wch1[wch1[, .I[which.medX(mean_int)], by=c("Protein","Feature","Mixture")]$V1]}
      
        wch2 <- wch[wch.frc, on=c("Protein","Feature","Mixture","Run")] 
        cnt.wch2 <- wch2[, c(.(cnt = sum(!is.na(Abundance)))), by = list(Mixture, Protein, Feature)] 
        cnt.wch0 <- wch0[, c(.(cnt = sum(!is.na(Abundance)))), by = list(Mixture, Protein, Feature)] 
          # table(cnt.wch2$cnt)
          # 1     2     3     4     5     6 
          # 16    68    75    88   190 19223 
          # table(cnt.wch0$cnt)
          # 0     1     2     3     4     5     6 
          # 275   247   504   584   894  2035 81507 
        work2 <- rbindlist(list(wch0[,-c("cnt")], wch2[,-c("cnt", "mean_int")] ))

   } else { work2 <- copy(work1) }
      
      
  ########################################################
  ################################### remove Features with only NAs in all channels
  cnt.work2 <- work2[, c(.(cnt = sum(!is.na(Abundance)))), by = list(Mixture, Protein, Feature)] 
    # table(cnt.work2$cnt)
    #   0      1      2      3      4      5      6 
    # 275    263    572    659    982   2225 100730 
  work3 <- work2[cnt.work2[cnt != 0,-"cnt"], on=c("Protein", "Mixture", "Feature")]

  cnt.work3 <- work3[, c(.(cnt = sum(!is.na(Abundance)))), by = list(Mixture, Protein, Feature)] 
  
    # tb <- table(cnt.work3$cnt)
    # dtb <- apply(data.matrix(as.data.frame(tb)) , 1, prod)
    # outperc <- round(100*(sum(dtb[1:5]) / dtb[6]),1)
    # print(paste0("By removing Features with not-complete Channel-sets, ",outperc, "% of the data will be taken out."))
    # if (!allChannelsIn) {"The filtering won't be done unless -allChannelsIn- set to TRUE."}
      
  ########################################################
  ################################### remove Features with not complete channel-set
  if (allChannelsIn) {
    cnt.work3 <- work3[, c(.(cnt = sum(!is.na(Abundance)))), by = list(Mixture, Protein, Feature)]
    
    work3 <- work3[cnt.work3[cnt == 6, -"cnt"], on=c("Protein", "Mixture", "Feature")]
    cnt.work3 <- work3[, c(.(cnt = sum(!is.na(Abundance)))), by = list(Mixture, Protein, Feature)] 
  }

    
  ########################################################
  ################################### continue with Proteins exist in all Mixrture    
    if (onlySharedPRT) {
      nPrtMix <- work3[, c(.(cnt = uniqueN(Mixture))), by = list(Protein)]
      work3 <- work3[Protein %in% nPrtMix[cnt==uniqueN(work3$Mixture), Protein] ]
      # work3 <- work3[nPrtMix[cnt==uniqueN(work3$Mixture), -"cnt"], on="Protein"]
      }
    
    
  ########################################################
  ################################### continue with Proteins with a mimimum number of Peptides in Mixture (minPptMIX)   
    if (minPptMIX > 1) {
      nPptMix <- work3[, c(.(cnt = uniqueN(Peptide))), by = list(Protein, Mixture)]
      work3 <- work3[Protein %in%  nPptMix[cnt>1, Protein]]
      # work3 <- work3[nPptMix[cnt>1, -"cnt"], on="Protein"]
    } 
    
    
  ################################################################# 
  ## VSN normalization  
    work.nrm <- copy(work3)
    work.nrm[, Intensity := 2^Abundance]
      if (toupper(FractComb) != "NONE") { # Fractions (Runs) are combined, continue with Mixtures
      work.nrm.w  <- data.table::dcast(work.nrm, Protein + Peptide + Charge ~ Mixture + Channel, value.var = "Intensity")
      cols.ch <- c("Protein","Peptide","Charge")
      cols <- names(work.nrm.w[,-cols.ch, with=FALSE])
      } else { # Fractions (Runs) are not combined
        work.nrm.w  <- data.table::dcast(work.nrm, Protein + Peptide + Charge + Run ~ Mixture + Channel, value.var = "Intensity")
        cols.ch <- c("Protein","Peptide","Charge", "Run")
        cols <- names(work.nrm.w[,-cols.ch, with=FALSE])     
        } 
        
        ## number of NAs per column and row
        na.row <- work.nrm.w[, Reduce(`+`, lapply(.SD, function(x) is.na(x))), .SDcols=cols] #NA per row    table(na.row)
        na.col <- work.nrm.w[, lapply(.SD, function(x) sum(is.na(x))), .SDcols = cols]       #NA per col    table(na.row)
        
        ## apply VSN 
        mat <- as.matrix(work.nrm.w[, cols, with=FALSE]) 
        fit <- vsn2(mat) 
        pred <- predict(fit, newdata = mat, useDataInFit = TRUE) 
        vsntest <-  meanSdPlot(pred, plot=FALSE)$gg
        
           as.data.table(pred) %>%
            cbind(work.nrm.w[,cols.ch, with=FALSE],.) %>%
             melt(., id.vars = cols.ch,
                 variable.name="factors",
                 value.name="Abundance.norm") %>%
                .[, c("Mixture", "Channel") := tstrsplit(factors, "_", fixed=TRUE)] %>%
                 .[,factors:=NULL] %>%
                  .[work.nrm, on=names(.)[!names(.) %in% "Abundance.norm"]] %>%
                   .[,Intensity:=NULL] -> work
                     rm(work.nrm)
        
        setnames(work, old="Abundance", new="Abundance.noNorm")
        setnames(work, old="Abundance.norm", new="Abundance")
        # setnames(work.NORM, old="Protein", new="Protein")
        # work.NORM$Mixture <- "Single"
        changeCols<- c(names(Filter(is.character, work)), names(Filter(is.integer, work)))
        work[,(changeCols):=lapply(.SD, as.factor),.SDcols=changeCols]
        
    
    # ## with strata
    #   work.nrm.w2 <- data.table::dcast(work.nrm, Protein + Peptide + Charge + Run + Mixture ~ Channel, value.var = "Intensity")
    #   cols2 <- names(work.nrm.w2[,-c("Protein","Peptide","Charge","Run")])
    #     
    #     mat2 <- work.nrm.w2[, cols2, with=FALSE]
    #     mat2$Mixture <- as.integer(mat2$Mixture)
    #     mat2 <- as.matrix(mat2)
    #       fit2 <- vsn2(mat2[,-1] , strata=as.integer(mat2[,1])) 
    #       pred2 <- predict(fit2, newdata = mat2[,-1], useDataInFit = TRUE)   
    #       fit3 <- vsn2(mat2[,-1]) 
    #       pred3 <- predict(fit3, newdata = mat2[,-1], useDataInFit = TRUE) 
    # 
    #       # plot(pred2,
    #       #      pred3,
    #       #      pch = ".", asp = 1, 
    #       #      col = hsv(seq(0, 1, length=36),
    #       #                0.8, 0.6)[mat2$Mixture],
    #       #      xlab = "without strata", 
    #       #      ylab = "print-tip strata")
    
    ## end VSN normalization  
    ################################################################# 

        
        
    ################################################################# 
    ## Median Polish 
    if (medpolON == "allMix") { # median polish over all mixture at once
      work.w <- data.table::dcast(work, Protein + Peptide + Charge + Run ~ Mixture + Channel, value.var = "Abundance")
      cols <- names(work.w)[(!colnames(work.w) %in% c("Protein", "Peptide", "Charge", "Run"))]
      work.mp  <- work.w[ , f.medpol(.SD, para=MPpara)[[1]] , by = .(Protein), .SDcols= cols  ]   # MPpara = "coleff" or "medres"
      work.mp[, c("Mixture", "Channel") := tstrsplit(Channel, "_",  fixed=TRUE)]
      # work.mp$MPmethod <- "allMix"
    }

    if (medpolON == "indMix"){ # median polish for each mixture separetly
      work.w <- data.table::dcast(work, Protein + Peptide + Charge + Mixture ~ Channel, value.var = "Abundance")
      cols <- names(work.w)[(!colnames(work.w.ind) %in% c("Protein", "Peptide", "Charge", "Mixture"))]
      work.mp <- work.w[ , f.medpol(.SD, para=MPpara)[[1]] , by = .(Protein, Mixture), .SDcols= cols  ]   # MPpara = "coleff" or "medres"
      # work.mp$MPmethod <- "indMix"
    }
    

    # work.mp <- rbind(work.mp.all, work.mp.ind)
    work.final <- work.mp[
      unique(work[, list(Protein, Mixture, Channel, Condition, BioReplicate)]),
      on=c("Protein", "Mixture", "Channel")]
    ################################################################# 
    
    
        ## break peptide names (see function .spltpep)
        tmp <- data.table(Peptide = as.character(unique(work$Peptide)), PeptideN=NA)
        tmp[, PeptideN := sapply(Peptide, function(x) .spltpep(x,n=9))]
        work <- work[tmp, on="Peptide"]
    
    
    ################################################################# 
    ## extract bioloical parameters  ..> genes
    
    # ens <- useMart("ensembl", "hsapiens_gene_ensembl")
    trans = getBM(attributes = c("uniprotswissprot", "mgi_symbol", "ensembl_gene_id", #"ucsc", "mgi_symbol"
                                 "external_gene_name",  "entrezgene", "description"), 
                  # filters = "uniprotswissprot", values = "Q64213", 
                  filters = "uniprotswissprot", values = unique(work.final$Protein),
                  uniqueRows=FALSE, mart = ens)
    trans <- as.data.table(trans)
    
    setkey(work.final, "Protein")
    setkey(trans, "uniprotswissprot")
    work.final[trans, Gene := external_gene_name]
    # work.final[Protein %in% "Q8CDM1" , Gene := "Atad2"]
    ################################################################# 
    
          ## test median polish -- plots
              # dd <- work[Protein == "Q6ZQK0"]
              # ddd <- work.final[Protein == "Q6ZQK0"]
              # plot.profile(dd,ddd)
              # dd <- work[Protein == "Q9QZ67"]
              # ddd <- work.final[Protein == "Q9QZ67"]
              # plot.profile(dd,ddd)             
              
              #
              # dd <- work[Protein=="Q9D8V7"]
              # ddd <- work.final[Protein=="Q9D8V7"]
              # 
              # dd <- work[Protein=="Q6P2L6"]
              # ddd <- work.final[Protein=="Q6P2L6"]
              # 
              # dd <- work[Protein=="Q64287"]
              # ddd <- work.final[Protein=="Q64287"]    
              # plot.profile(dd,ddd)
              # 
              # ddt <- work[Protein=="Q64287"]
              # prt.plot(ddt)
  
    # dd <- work[Protein == "Q3TTA7"]
    # ddd <- work.final[Protein == "Q3TTA7"]
    # plot.profile(dd,ddd)
              
 
    if (mix3rev) {
      
      work.final[Mixture ==2 & Condition == "WT.0h", Condition := "__Cblb.0h"]
      work.final[Mixture ==2 & Condition == "WT.24h", Condition := "__Cblb.24h"]
      work.final[Mixture ==2 & Condition == "WT.48h", Condition := "__Cblb.48h"]
      
      work.final[Mixture ==2 & Condition == "Cblb.0h", Condition := "WT.0h"]
      work.final[Mixture ==2 & Condition == "Cblb.24h", Condition := "WT.24h"]
      work.final[Mixture ==2 & Condition == "Cblb.48h", Condition := "WT.48h"]
      
      work.final$Condition <- gsub("__", "" ,work.final$Condition) 
      work.final$Condition <- as.factor(work.final$Condition)
      
    }
    
  ## this will be need  
  workf <- copy(work.final)
    
  ################################################################# 
  ## statistical analysis -- limma
  work.final[, BioReplicate:= do.call(paste, c(.SD, sep = "_")), .SDcols=c("Mixture", "Channel", "BioReplicate")]
  work.final$Run <- "single"
  work.final <-  unique(work.final)  
  
  
  ## define Contrast
    dt <- work.final
    dt$Protein <-  paste0(dt$Gene,"_", dt$Protein )
       ## dgenerate matrix
       dtw <- data.table::dcast(dt, Protein ~ BioReplicate, value.var = "Abundance")
       mat <- as.matrix(dtw[,-c("Protein")])
       rownames(mat) <- dtw[,Protein]
       
       reftb <- unique(dt[,c("BioReplicate", "Condition", "Run", "Channel")])
       reftb[, c("Mixture", "Channel", "Replicate") := tstrsplit(BioReplicate, "_",  fixed=TRUE)]

         Condition <- as.character(reftb$Condition)
         Mixture <- as.character(reftb$Mixture)
         Channel <- as.character(reftb$Channel)
         Replicate <- as.character(reftb$Replicate)
           
         design <- model.matrix(~ 0 + Condition)
         fit <- lmFit(mat, design)
         
         
         ########################### for pairwise comparison
             ## store inference results
             fitList <- NULL
             ## store top rankrd proteins
             topList <- NULL

             conds <- paste("Condition", unique(Condition), sep = "")
               for(j in 1:(length(conds)-1)){
                 for(k in (j+1):length(conds)){
                   cond1 <- conds[j]
                   cond2 <- conds[k]
                   comp <- paste(cond1, cond2, sep="-")
                   contrast <- makeContrasts(contrasts=comp, levels=design)
                  
                   fit2 <- contrasts.fit(fit, contrast)
                   fit2 <- eBayes(fit2)
                   
                     ## topTABLE to file
                       tab <- topTable(fit2, coef=1, adjust="BH", number=round((TopNperc/100)*nrow(mat)))
                       tab$Protein <- rownames(tab)
                       setDT(tab)
                       tab[, Comparison := gsub("Condition", "", comp)]
                       tab[, id := seq(nrow(tab))]
                       tab[, c("Gene", "Protein") := tstrsplit(Protein, "_",  fixed=TRUE)]
                       setcolorder(tab, c("id","Comparison","Protein","Gene","logFC","P.Value","adj.P.Val","AveExpr","t","B"))
                       topList <- rbind(topList, tab)
                     ##
                   
                       
                  fit.dt <- data.table(
                    id = seq(nrow(fit2$coefficients)),
                    Comparison = gsub("Condition", "", paste(conds[j], conds[k],sep="-")),
                    Protein = rownames(fit2$coefficients),
                    P.Value = as.vector(fit2$p.value),
                    logFC = as.vector(fit2$coefficients),
                    DF = as.vector(fit2$df.total),
                    SE = as.vector(sqrt(fit2$s2.post) * fit2$stdev.unscaled)
                  )
                  fit.dt[, c("Gene", "Protein") := tstrsplit(Protein, "_",  fixed=TRUE)]
                  fit.dt$adj.P.Val <- p.adjust(fit.dt$P.Value, method="BH")
                  
                  fitList <- rbind(fitList, fit.dt)
                 }
               }
             
             # fitList$adj.P.Val.allPairs <- p.adjust(fitList$P.Value, method="BH")
             setcolorder(fitList, c("id","Comparison","Protein","Gene","logFC","P.Value","adj.P.Val","DF","SE"))
         ########################### 
            
             
             
             
         ########################### for SINGLE comparison (all WT v. all KO)
         contrast <- t(matrix(c(1,1,1,-1,-1,-1), nrow=1))
             # # model by Condition + Mixture
             # contrast <- t(matrix(c(1,1,1,-1,-1,-1,0,0), nrow=1))
         rownames(contrast) <- colnames(fit$coefficients)
         colnames(contrast) <- comp <- "WT-CBLB"
         
           fit2 <- contrasts.fit(fit, contrast)
           fit2 <- eBayes(fit2)
           
             ## topTABLE to file
             tab.full <- topTable(fit2, coef=1, adjust="BH", number=round((TopNperc/100)*nrow(mat)))
             tab.full$Protein <- rownames(tab.full)
             setDT(tab.full)
             tab.full[, Comparison := gsub("Condition", "", comp)]
             tab.full[, id := seq(nrow(tab.full))]
             tab.full[, c("Gene", "Protein") := tstrsplit(Protein, "_",  fixed=TRUE)]
             setcolorder(tab.full, c("id","Comparison","Protein","Gene","logFC","P.Value","adj.P.Val","AveExpr","t","B"))
             topList <- rbind(tab.full, topList)
             ##
           
           
           fit.dt.full <- data.table(
             id = seq(nrow(fit2$coefficients)),
             Comparison = comp,
             Protein = rownames(fit2$coefficients),
             P.Value = as.vector(fit2$p.value),
             logFC = as.vector(fit2$coefficients),
             DF = as.vector(fit2$df.total),
             SE = as.vector(sqrt(fit2$s2.post) * fit2$stdev.unscaled)
           )
           fit.dt.full[, c("Gene", "Protein") := tstrsplit(Protein, "_",  fixed=TRUE)]
           fit.dt.full$adj.P.Val <- p.adjust(fit.dt$P.Value, method="BH")
           setcolorder(fit.dt.full, c("id","Comparison","Protein","Gene","logFC","P.Value","adj.P.Val","DF","SE"))
           fitList <- rbind(fit.dt.full, fitList)
           ########################### 
           

     ################################################################# 
     ## Generate PCA plots
           
     prtvar <- workf[, c(.(var = var(Abundance))), by = list(Protein)]
     tmp <- prtvar[order(-var)]
     tmp <- tmp[1:500,]
     tmp <- workf[Protein %in% tmp$Protein,]

     wk <- data.table::dcast(tmp, Protein ~ Condition + Mixture + BioReplicate, value.var = "Abundance")
     wk$na_count <- apply(wk[,-1], 1, function(x) sum(is.na(x)))
     wk <- wk[na_count ==0]
     mat <- as.matrix(wk[,-c("Protein","na_count")])
     pca2 = prcomp(t(mat), scale = TRUE)

     pcaplot2 <-
       fviz_pca_ind(pca2,
                    col.ind = "cos2", # Color by the quality of representation
                    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                    repel = TRUE)

          
     
    # tmp <- fitList[, head(.SD, 500), by=list(Comparison), .SDcols= "SE"] 
    # tmp[, idx := 1:.N, by = Comparison]
    # tmpw <- data.table::dcast(tmp, idx~Comparison, value.var="SE")
    #   tmpw$na_count <- apply(tmpw[,-1], 1, function(x) sum(is.na(x)))
    #   tmpw <- tmpw[na_count ==0]
    #   mat1 <- as.matrix(tmpw[,-c("idx", "na_count")])
    #   pca11 = prcomp(t(mat1), scale = TRUE)           
    #   
    #   pcaplot1 <- 
    #     fviz_pca_ind(pca11, 
    #                  col.ind = "cos2", # Color by the quality of representation
    #                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
    #                  repel = TRUE)     # Avoid text overlapping           
    #        
           
     wk <- data.table::dcast(workf, Protein ~ Condition + Mixture + BioReplicate, value.var = "Abundance")
     wk$na_count <- apply(wk[,-1], 1, function(x) sum(is.na(x)))
     wk <- wk[na_count ==0]
     mat <- as.matrix(wk[,-c("Protein","na_count")])
     pca1 = prcomp(t(mat), scale = TRUE)
     # pca2 = princomp(mat, cor = TRUE)
     library(factoextra)
     
     pcaplot <- 
       fviz_pca_ind(pca1, 
                    col.ind = "cos2", # Color by the quality of representation
                    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                    repel = TRUE)     # Avoid text overlapping
           
           
           
           
           
           
     ################################################################# 
     ## Generate Heatmaps            
     
     ## get the top N Proteins of each comparison
     topN <- 10
     topPrtCom = topList[topList[, .I[1:topN], Comparison]$V1, list(Comparison, Protein, Gene)]   
       hm.protein <- unique(topPrtCom$Protein)
     
       
     hm.dt <- work.final[Protein %chin% hm.protein] 
       if (!"Cblb" %in% hm.dt$Gene) { 
       hm.dt <- rbind(work.final[Gene %chin% c("Irf4","Cblb")], hm.dt) ## add cblb and Irf4
       }
     
     hm.dt[, Rep := tstrsplit(BioReplicate, "_",  fixed=TRUE)[3]]
     hm.dt[, xTag := do.call(paste, c(.SD, sep = "\nRep")), .SDcols=c("Condition", "Rep")]
     hm.dt[, yTag := do.call(paste, c(.SD, sep = "_")), .SDcols=c("Protein", "Gene")]
     
     hm.dt.w <- data.table::dcast(hm.dt, yTag ~ xTag, value.var="Abundance" )
     hm.mat <- as.matrix(hm.dt.w[,-1])
     rownames(hm.mat) <- hm.dt.w$yTag
     
# 
#      mypal <- colorRampPalette( brewer.pal( 10 , "RdYlBu" ) )
#      hm.plot <- ggplot(data=hm.dt, aes(x=xTag, y=yTag, fill=Abundance)) +
#        geom_tile(color="white", size=0.05) +
#        scale_fill_distiller(palette = "Spectral",na.value = "grey70") +
#        facet_grid(~Mixture, scales = "free", space = "free", switch="y") +
#        scale_fill_manual( values = mypal )+
#        theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
     
     
     library("RColorBrewer")
     
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
     
     hmallclust.plot<- Heatmap(hm.mat, name = "Abundance", 
                               # column_dend_height = unit(2, "cm"), 
                               row_dend_width = unit(2, "cm"),
                               na_col = "white", 
                               clustering_distance_rows = "pearson",
                               cluster_columns = TRUE,
                               # col = circlize::colorRamp2(c(3,7), c("darkblue", "yellow"))
                               col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256),
                               top_annotation = ha)
     
     
                 ## save toptable and fitlists
                if (mix3rev) {
                  savePhrase <- paste0( "__mix3rev_",
                                        "__FracComb_",FractComb,
                                        "__rmNAchnnl_",allChannelsIn,"__whichPSM_",whichPSM,
                                        "__MedPol_",medpolON,"__measurements_",meas,
                                        "__rmMix2_",remMix2, "__onlySharedPRT_",onlySharedPRT,
                                        "__minPptMIX_",minPptMIX)
                } else {
                 savePhrase <- paste0( "__FracComb_",FractComb,
                                       "__rmNAchnnl_",allChannelsIn,"__whichPSM_",whichPSM,
                                       "__MedPol_",medpolON,"__measurements_",meas,
                                       "__rmMix2_",remMix2, "__onlySharedPRT_",onlySharedPRT,
                                       "__minPptMIX_",minPptMIX)
                }
     
                 write.csv(topList, file = paste0(getwd(),"/loop_all/tophits_data/topTable_limma", savePhrase, ".csv"))
                 write.csv(fitList, file = paste0(getwd(),"/loop_all/statistics_data/fitList_limma", savePhrase, ".csv"))
                 
                 savefile <- paste0(getwd(),"/loop_all/volcano/PairwiseVolcanoPlots_",savePhrase,".pdf")
                 volc_cblb(fitList, savefile=savefile)   
                 savefile <- paste0(getwd(),"/loop_all/volcano_notext/PairwiseVolcanoPlots_",savePhrase,".pdf")
                 volc_cblb(fitList, savefile=savefile, label = TRUE)    
                 
                 
                 pdf(paste0(getwd(),"/loop_all/heatmap/TopProt_Heatmap_",savePhrase,".pdf"), width = 12, height = 11)
                 print(hm.plot)
                 dev.off()
                 
                 pdf(paste0(getwd(),"/loop_all/heatmap_allclust/TopProt_Heatmap_",savePhrase,".pdf"), width = 12, height = 11)
                 print(hmallclust.plot)
                 dev.off()
                 
                 pdf(paste0(getwd(),"/loop_all/pca/PCA_",savePhrase,".pdf"), width = 9, height = 7)
                 print(pcaplot)
                 dev.off()                 
           
              }
            }
          }
        }
      }
    }
#   }
# }

    



 
 # # comparison<-matrix(c(-1,-1,-1,1,1,1), nrow=1)
 # workf$Protein <-  paste0(workf$Gene,"_", workf$Protein )
 # limma.sig_newMeas <- MSstatsTMT::groupComparison.TMT(as.data.frame(workf),
 #                                                      contrast.matrix = 'pairwise',
 #                                                      # contrast.matrix = comparison,
 #                                                      remove_norm_channel = FALSE,
 #                                                      model = 'limma',
 #                                                      adj.method = "BH")
 # volc_cblb(limma.sig_newMeas, main="limma_all_noFilt_v2_cblbVwt_allFrac_noMix2")
 # write.csv(limma.sig_newMeas, file = "limma_sig_noFilt_v1.csv")
 

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ###########################################
  ###########################################
  ###########################################
  # @Andreas: This block is supposed to be taken out. This is the old workflow. The outcome might be different from how it
  # supposed to be, as the data format has almost completely changed in the last few days. I will use some parts for the
  # new algorithm, but it needs to be recoded to output desired results.
  
    ## decision1 : use rows with most number of non-NA measurement
    cnt.tmp1 <- work.tmp[, c(.(cnt = sum(!is.na(Intensity)))), by = list(Run, Feature, Channel)]
    cnt.tmp1 <- cnt.tmp1[cnt.tmp1[, .I[cnt == max(cnt, rm.na=TRUE)], by=c("Run","Feature", "Channel", "Ion.Score")]$V1]
    if (sum(cnt.tmp1$cnt > 1)) {
      work.tmp1 <- work.tmp[cnt.tmp1, on=c("Run", "Feature", "Channel", "Ions.Score")] }
    
      ## decision2 : keep the row with higher identification score
      cnt.tmp2 <- work.tmp1[, c(.(cnt = uniqueN(Ions.Score))), by = list(Run, Feature, Channel)] 
      if (sum(cnt.tmp2$cnt > 1)) {
      work.tmp2 <- work.tmp1[work.tmp1[, .I[Ions.Score == max(Ions.Score, rm.na=TRUE)], by=c("Run","Feature", "Channel")]$V1]}
      
        ## decision2 : keep the row with max intensity
        cnt.tmp3 <- work.tmp2[, c(.(cnt = length(Intensity))), by = list(Run, Feature, Channel, Ions.Score)]
        if (sum(cnt.tmp3$cnt > 1)) {
        work.tmp3 <- work.tmp2[work.tmp2[, .I[which.max(Intensity)], by=c("Run","Feature", "Channel")]$V1]}
        
          ## merging data        
          work.tmp4 <- rbindlist(list(work.unq, work.tmp3[,-c("i.cnt","cnt")] ))
          work.tmp4[, PSM := do.call(paste, c(.SD, sep = "_")), .SDcols=c("Peptide", "Charge")]
          work.tmp4[, Feature := NULL] 
          
            ## remove features (PSM) which has missing measurements within each run (less than number of chennels)
            missPSM <- work.tmp4[, c(.(cnt = sum(!is.na(Intensity)))), by = list(Run, Protein, PSM)] [cnt==uniqueN(work$Channel),]
            work.tmp5 <- work.tmp4[missPSM, on=c("Protein", "Run", "PSM")] 
            
              ## Remove the peptide ions (PSM) overlapped among multiple fractions (Run) of same biological mixture (Mixture)
              # sharedPSM.Mix <- work.tmp5[, c(.(cnt = uniqueN(Run))), by = list(PSM, Mixture)][cnt==1,]
              # work.tmp6 <- work.tmp5[sharedPSM.Mix, on=c("PSM", "Mixture")] 
            
              ##!!! UPDATE 14.06.2015: Marc suggested that instead of removing overlapped PSM between Runs of
              ##!!! the same biological mixture, take either Maximum or Sum the redundant values.
              sharedPSM.Mix <- work.tmp5[, c(.(cnt = uniqueN(Run))), by = list(PSM, Mixture)]
                wch0 <- work.tmp5[sharedPSM.Mix[cnt==1,], on=c("PSM", "Mixture")] 
                wch  <- work.tmp5[sharedPSM.Mix[cnt>1,], on=c("PSM", "Mixture")] 
                  wch1 <- wch[ , .(mean_int = mean(Intensity, na.rm=TRUE)), by = c("Protein","PSM","Mixture","Run")]
                  wch1 <- wch1[wch1[, .I[which.max(mean_int)], by=c("Protein","PSM","Mixture")]$V1]
                    wch2 <- wch[wch1, on=c("Protein","PSM","Mixture","Run")] 
                      cnt.wch2 <- wch2[, c(.(cnt = sum(!is.na(Intensity)))), by = list(Mixture, Protein, PSM)] [cnt==uniqueN(work$Channel),]
                      wch3 <- wch2[cnt.wch2, on=c("Protein", "Mixture", "PSM")] 
                      work.tmp6 <- rbindlist(list(wch0[,-c("cnt")], wch3[,-c("cnt", "mean_int", "i.cnt")] ))
                      
              # test1 <- work.tmp6[, c(.(cnt = uniqueN(Run))), by = list(PSM, Mixture)]
              # test2 <- wch.f[, c(.(cnt = uniqueN(Run))), by = list(PSM, Mixture)]
              # table(test1$cnt)
              # table(test2$cnt)
              # #@@ these are for checking
              # wch2[PSM %in% "eTVSEESNVLcLSk_3" & Mixture==2 ,]
              # a <- cnt.tmp1[Run == 39 & Feature %in% "P08249_vDFPQDQLATLTGR_2",]
              # a.cnt <- a[, c(.(cnt = uniqueN(Channel))), by = list(Ions.Score)]


                ## remove single-shot proteins (here I count PSM, counting peptide might be better!!)
                snglProt <- work.tmp6[, c(.(cnt = uniqueN(Peptide))), by = list(Protein)][cnt!=1, Protein]
                work.tmp7 <- work.tmp6[Protein %in% snglProt,] ## if unique PEP: 2536 PRT -- if unique PSM: 2519 PRT
              
                
                
## dataset clean ---> after applying filters
work.clean <- work.tmp7[!is.na(Protein),]

## Each Run is corresponding to a Fraction. We combine the fractions of each biological mixture with each other and 
## attribute the mixture id (2,3,4) to Run. 
work.clean[, Fraction.id := Run]
work.clean[, Run := Mixture]

work.clean[, Rep:= do.call(paste, c("Rep",.SD,sep = ".")), .SDcols=c("BioReplicate")]
work.clean[, Mix:= do.call(paste, c("Mix",.SD,sep = ".")), .SDcols=c("Mixture")]
work.clean[, BioReplicate:= do.call(paste, c(.SD, sep = ".")), .SDcols=c("Run","Channel","BioReplicate")]

cols <- c("Mix","Rep","BioReplicate")
work.clean[, (cols):=lapply(.SD, factor), .SDcols=cols]

## log2 transformation
work.clean$log2Intensity <- log2(work.clean$Intensity) 

## Number of negative values : if intensity is less than 1, replace with zero
## then we don't need to worry about -Inf = log2(0) !!!!!!!!!!!!! do it.
###########################################
###########################################
###########################################



###################################### Summarization to prtein level with MedianPolish
## Function to apply median polish to data table blocks. This is much faster than for loops.

work.w <- data.table::dcast(work.clean, Protein + PSM + Run ~ Channel, value.var = "log2Intensity")[,-c("PSM")]
cols <- names(work.w)[(!colnames(work.w) %in% c("Protein", "Run"))]
start_time <- Sys.time()
  work.sum <- work.w[ , f.medpol(.SD) , by = .(Protein, Run), .SDcols= cols  ]   # simplest
end_time <- Sys.time()
end_time - start_time


## Merging summarized data with clean dataset
work.clean.sub <- work.clean[, c("Run","Protein", "Channel", "BioReplicate", "Condition",  "Mix", "Rep")]
work.clean.sub <-  unique(work.clean.sub)
work.preNORM <- work.sum[work.clean.sub, on=c("Protein", "Run", "Channel")]

## check
# prt <- c("Q9CY57", "P17710", "Q6ZPJ3")
# a <- work.clean[Protein %in% prt,]
# aw <- data.table::dcast(a, Protein + PSM + Run ~ Channel, value.var = "log2Intensity")[,-c("PSM")]
# cols <- names(aw)[(!colnames(aw) %in% c("Protein", "Run"))]
# at <- aw[ , f.medpol(.SD) , by = .(Protein, Run), .SDcols= cols  ]   # simplest
# atw <- data.table::dcast(at[Protein %in% "P17710",], Protein + Run ~ Channel, value.var = "Abundance")
# qnt <- data.table::dcast(quant.byprotein.NN[Protein %in% "P17710",], Protein + Run ~ Channel, value.var = "Abundance")
# cblb protein  Q13191


# dt1 <- work.preNORM[Mixture == 2]
# dt2 <- work.preNORM[Mixture == 3]
# dt3 <- work.preNORM[Mixture == 4]
# prtlist <- Reduce(intersect, list(dt1$Protein, dt2$Protein, dt3$Protein))
# work.preNORM <- work.preNORM[Protein %chin% prtlist]
# 
# 
# (a <- structable(Condition  ~ Channel, data = psm.ano))






# ## Significance inference
# data <- as.data.frame(work.NORM)
# colnames(data)[colnames(data) == 'BioReplicate'] <- 'Subject'
# colnames(data)[colnames(data) == 'Condition'] <- 'Group'
# groups <- unique(work.NORM$Condition)
# 
# comparison<-matrix(c(-1,-1,-1,1,1,1),nrow=1)
# 
# row.names(comparison)<- c("groupCblb.0h-groupCblb.24h-groupCblb.48h", "groupWT.0h-groupWT.24h-groupWT.48h" )
# colnames(comparison)<- groups
# 
# limma.sig_new <- groupComparison.TMT(as.data.frame(work.NORM),
#                                  contrast.matrix = 'pairwise',
#                                  # contrast.matrix = comparison,
#                                  remove_norm_channel = FALSE,
#                                  model = 'limma',
#                                  adj.method = "BH")
# volc_cblb(limma.sig_new, main="limma_all_new")
# write.csv(limma.sig_new, file = "limma_sig_new.csv")





limma.sig_newMeas <- groupComparison.TMT(as.data.frame(work.NORM),
                                        contrast.matrix = 'pairwise',
                                        # contrast.matrix = comparison,
                                        remove_norm_channel = FALSE,
                                        model = 'limma',
                                        adj.method = "BH")
volc_cblb(limma.sig_newMeas, main="limma_all_noFilt_v1")
write.csv(limma.sig_newMeas, file = "limma_sig_noFilt_v1.csv")






## plots
prt <- "A0A023T778"
prt <- "A0A075B5P2"


inp <- work[Protein == prt,]
qnt <- work.preNORM[Protein %chin% prt,]
inp[,Run:=NULL]

setnames(inp, old = "Protein", new = "Protein")
setnames(inp, old = "Mixture", new = "Run")
inp[,id:=seq_along(Channel), by=c("Run")]

inp <- setorderv(inp, c("id", "Channel", "Run"), c(1,1,1))
qnt <- setorderv(qnt, c("Run", "Channel"), c(1,1))


ggplot(inp, aes(x=id, y=log2(Intensity), col=Peptide)) + 
  geom_point(aes(shape=as.factor(Charge))) +  facet_grid(Run~Condition) + 
  geom_line(aes(linetype=as.factor(Charge)), stat="identity") +
  theme(axis.text.x = element_text(size=8, angle = 90, vjust = 0.5))

ggplot(work, aes(x=))
