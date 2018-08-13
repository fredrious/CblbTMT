

################################################################# 
## statistical analysis -- limma
## returns parirwise and contrast defind statistics

statSig <- function(workf, TopNperc=5) {
  
  dt <- copy(workf)
  dt[, BioReplicate:= do.call(paste, c(.SD, sep = "_")), .SDcols=c("Mixture", "Channel", "BioReplicate")]
  dt$Run <- "single"
  dt <-  unique(dt)  
  
  
  ## define Contrast
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
  
  conds <- sort(factor(paste("Condition", unique(Condition), sep = "")))
  
  for(j in 1:(length(conds)-1)){
    for(k in (j+1):length(conds)){
      cond1 <- conds[j]
      cond2 <- conds[k]
      comp <- factor(paste(cond1, cond2, sep="-"))
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
        Comparison = factor(gsub("Condition", "", paste(conds[j], conds[k],sep="-"))),
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
  colnames(contrast) <- comp <- factor("WT-Cblb")
  
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
  fit.dt.full$adj.P.Val <- p.adjust(fit.dt.full$P.Value, method="BH")
  setcolorder(fit.dt.full, c("id","Comparison","Protein","Gene","logFC","P.Value","adj.P.Val","DF","SE"))
  fitList <- rbind(fit.dt.full, fitList)
  
  return(list(topList, fitList))
}

