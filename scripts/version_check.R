


view.volc <- function (tag) {
  loadPhrase <- as.character(key.dt[keytag %in% tag, file])
  file <- paste0("/Users/farhad/_Rspace/_prj_Wolf_Isabelle/analysis/loopAll/volcano/volc_all_",loadPhrase,".pdf")
  system(paste0('open "', file, '"'))
}


view.hmWTKO <- function (tag) {
  loadPhrase <- as.character(key.dt[keytag %in% tag, file])
  file <- paste0("/Users/farhad/_Rspace/_prj_Wolf_Isabelle/analysis/loopAll/heatmap/hmapWTKO_",loadPhrase,".pdf")
  system(paste0('open "', file, '"'))
}


view.hmAll <- function (tag) {
  loadPhrase <- as.character(key.dt[keytag %in% tag, file])
  file <- paste0("/Users/farhad/_Rspace/_prj_Wolf_Isabelle/analysis/loopAll/heatmap/hmapALL_",loadPhrase,".pdf")
  system(paste0('open "', file, '"'))
}


Frac10.arr = c(TRUE, FALSE) 
FractComb.arr <-  c("max", "single")    
allChannelsIn.arr <- c(TRUE, FALSE)     
remMix2.arr <- c(FALSE, TRUE)           
minPptMIX.arr <- c(0,1,2) 
betweenFracNorm.arr <- c("median", "quantile", "none")
betweenChnNorm.arr = c("median","vsn", "qntVSN.noCal", "medVSN.noCal")

key.dt <- data.table(NULL)

for (ix in 1:length(Frac10.arr)) {
  for (jx in 1:length(FractComb.arr)) {
    for (kx in 1:length(allChannelsIn.arr)) {
      for (lx in 1:length(remMix2.arr)) {
        # for (mx in 1:length(onlySharedPRT.arr)) {
        for (nx in 1:length(minPptMIX.arr)) {
          for (px in 1:length(betweenFracNorm.arr)) {
            for (qx in 1:length(betweenChnNorm.arr)) {
              
            keytag <- paste(c(ix, jx, kx, lx, nx, px, qx), collapse = ".") 
              
            Frac10          <- Frac10.arr[ix]
            FractComb       <- FractComb.arr[jx]
            allChannelsIn   <- allChannelsIn.arr[kx]
            remMix2         <- remMix2.arr[lx]
            minPptMIX       <- minPptMIX.arr[nx]
            betweenFracNorm <- betweenFracNorm.arr[px]
            betweenChnNorm  <- betweenChnNorm.arr[qx]
            
            savePhrase <- paste0( keytag,"__Frac10_",Frac10,"__Frac_",FractComb,
                                  "__rmNAch_",allChannelsIn,
                                  "__rmM2_",remMix2,
                                  "__minPpt_",minPptMIX, "__FracNrm_", betweenFracNorm, 
                                  "__ChnNrm_", betweenChnNorm) 
            
            savePhrase <- gsub("TRUE","T", savePhrase)
            savePhrase <- gsub("FALSE","F", savePhrase)
            
            
            key.tmp <- data.table(keytag = keytag,
                                  Frac10 = Frac10[ix],
                                  FractComb = FractComb.arr[jx],
                                  allChannelsIn = allChannelsIn.arr[kx],
                                  remMix2 = remMix2.arr[lx],
                                  minPptMIX = minPptMIX.arr[nx],
                                  betweenFracNorm = betweenFracNorm.arr[px],
                                  betweenChnNorm = betweenChnNorm.arr[qx],
                                  file = savePhrase
            )

            key.dt <- rbind(key.dt, key.tmp)
              
                
            }
          }
        }
      }
    }
  }
} 

key.dt[, paste0("v",1:7) := tstrsplit(keytag, ".",  fixed=TRUE)]
key.dt[v1 == 2, Frac10 := FALSE ]
key.dt <- key.dt[,-c("v1","v2","v3","v4","v5","v6","v7")]
       
       
dir <- "/Users/farhad/_Rspace/_prj_Wolf_Isabelle/analysis/loopAll/limmaTops/"
topfl <- list.files(path = dir, full.names=FALSE)
iv <- gsub("Top_", "", topfl)
iv <- gsub("\\__.*", "", iv, fixed=FALSE)

dtops <- lapply(paste0(dir, topfl), read.csv)
names(dtops) <- iv

dtop <- rbindlist(dtops, fill=TRUE, idcol=TRUE)
dtop <- char2fact(dtop)
setnames(dtop, old=".id", new="keytag")



## criterion 1: P.Value of Cblb
# subtop <- dtop[Gene=="Cblb" & Comparison %like% "Cblb" & Comparison %like% "WT",
#                c("keytag", "Comparison", "P.Value","logFC", "adj.P.Val")]
subtop <- dtop[Gene %in% c("Cblb") & Comparison == "Cblb.24h-WT.24h",
               c("keytag", "Comparison", "P.Value","logFC", "adj.P.Val")]

sub <- merge(subtop, key.dt[,-c("file")], by="keytag", all.x = TRUE)
sub <- sub[minPptMIX == 1 & !betweenFracNorm %in% "quantile" & remMix2 == FALSE,]
sub <- sub[order(adj.P.Val)]
sub[1:100,]


## criterion 2: count Sig Genes in WT-Cblb
# sigdt <- dtop[adj.P.Val < 0.05 & abs(logFC) > 1 & Comparison %like% "WT" & Comparison %like% "Cblb"]
sigdt <- dtop[adj.P.Val < 0.05 & abs(logFC) > 1 & Comparison %like% "Cblb.24h-WT.24h" ]
sigdt <- sigdt[, .(cntSig = length(Protein)), by=keytag]
sigdt <- merge(sigdt, key.dt[,-c("file")], by="keytag", all.x = TRUE)
sigdt <- sigdt[remMix2 == FALSE & minPptMIX > 0 & !betweenFracNorm %in% "quantile"]
sigdt <- sigdt[order(-cntSig)]
sigdt[1:100,]

x <- sigdt[1:100,]
y <- sub[1:100,]
z <- merge(x, y, by="keytag", all = FALSE)
z


view.hmWTKO("2.2.2.1.2.1.3")
view.hmAll("2.2.2.1.2.1.3")
view.volc("2.2.2.1.3.3.2")



