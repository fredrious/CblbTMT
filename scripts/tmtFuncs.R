
################################
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
################################



###################################
## change Charachter and Integer columns to Factor
char2fact <- function(dt){
  changeCols<- c(names(Filter(is.character, dt)), names(Filter(is.integer, dt)))
  dt[,(changeCols):=lapply(.SD, as.factor),.SDcols=changeCols]
return(dt)
}

int2fact <- function(dt){
  changeCols<- c(names(Filter(is.character, dt)), names(Filter(is.integer, dt)))
  dt[,(changeCols):=lapply(.SD, as.factor),.SDcols=changeCols]
  return(dt)
}
###################################




################################
## function: median polish
f.medpol <-  function(dws, para) {
  clnm <- colnames(dws)
  dw <- as.matrix(dws)
  mp  <-  stats::medpolish(dws, na.rm=TRUE, trace.iter = FALSE, maxiter = 30)
  if (toupper(para)=="MEDRES") {
    tmp <- mp$overall + apply(mp$residuals, 2, function(x) mean(x, na.rm = TRUE)) # Abundance = overall median + mean(residual)
  } else { tmp <- mp$overall + mp$col } # Abundance = overall median + column effect
  result <- data.table(Channel=clnm, Abundance=tmp)
  res <- mp$residual
  return(list(result,res))
}


MedPolAll <- function(dw, MPpara="coleff") {
  dw.w <- data.table::dcast(dw, Protein + Peptide + Charge ~ Mixture + Channel, value.var = "Abundance")
  cols <- names(dw.w)[(!colnames(dw.w) %in% c("Protein", "Peptide", "Charge", "Run"))]
  dw.mp  <- dw.w[ , f.medpol(.SD, para=MPpara)[[1]] , by = .(Protein), .SDcols= cols  ]   # MPpara = "coleff" or "medres"
  dw.mp[, c("Mixture", "Channel") := tstrsplit(Channel, "_",  fixed=TRUE)]
  # work.mp$MPmethod <- "allMix"
  return(dw.mp)
}


MedPolInd <- function(dw, MPpara="coleff") {
  dw.w <- data.table::dcast(dw, Protein + Peptide + Charge + Mixture ~ Channel, value.var = "Abundance")
  cols <- names(dw.w)[(!colnames(dw.w) %in% c("Protein", "Peptide", "Charge", "Mixture"))]
  dw.mp <- dw.w[ , f.medpol(.SD, para=MPpara)[[1]] , by = .(Protein, Mixture), .SDcols= cols  ]   # MPpara = "coleff" or "medres"
  # work.mp$MPmethod <- "indMix"
  return(dw.mp)
}
################################




################################
## VSN normalization
VSNnorm <-  function(wdn, calib) {
  
  wdn[, Intensity := 2^Abundance]
  if (toupper(FractComb) != "NONE") { # Fractions (Runs) are combined, continue with Mixtures
    wdn.w  <- data.table::dcast(wdn, Protein + Peptide + Charge  ~ Mixture + Channel, value.var = "Intensity")
    cols.ch <- c("Protein","Peptide","Charge")
    cols <- names(wdn.w[,-cols.ch, with=FALSE])
  } else { # Fractions (Runs) are not combined
    wdn.w  <- data.table::dcast(wdn, Protein + Peptide + Charge + Run ~ Mixture + Channel, value.var = "Intensity")
    cols.ch <- c("Protein","Peptide","Charge", "Run")
    cols <- names(wdn.w[,-cols.ch, with=FALSE])     
  } 
  
  ## apply VSN 
  mat <- as.matrix(wdn.w[, cols, with=FALSE]) 
  fit <- vsn2(mat, calib=calib) 
  pred <- predict(fit, newdata = mat, useDataInFit = TRUE) 
  vsntest <-  meanSdPlot(pred, plot=FALSE)$gg + theme_bw() +  scale_fill_distiller(palette = "Spectral") 
  
  as.data.table(pred) %>%
    cbind(wdn.w[,cols.ch, with=FALSE],.) %>%
    melt(., id.vars = cols.ch,
         variable.name="factors",
         value.name="Abundance") %>%
    .[, c("Mixture", "Channel") := tstrsplit(factors, "_", fixed=TRUE)] %>%
    .[, factors := NULL] %>%
    .[wdn, on=names(.)[!names(.) %in% "Abundance"]] %>%
    .[, Intensity := NULL] -> work
  
  return(list(work, vsntest))
}
################################




################################
## Function to get the median value
which.medX <-  function(x) {
  xx <-  x >= median(x, na.rm = TRUE)
  which( x == x[xx][which.min(x[xx])] )}
## gives the median if length(x)=2*n+1 or the larger of the two middle values if length(x)=2*n. 
################################




################################
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
################################





