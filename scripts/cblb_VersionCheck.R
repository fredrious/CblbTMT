

meas.arr <- c("old" ,"new")                   ## Which measurement
PSM2PPT.arr = c(TRUE, FALSE)                  ## PSM to Peptide (Ratio > Median) 
FractComb.arr <-  c("med", "max")             ## combine fractions
remMix2.arr <- c(TRUE, FALSE)                 ## remove mixture 2 or not (high missingness in mix2)
allChannelsIn.arr <- c(TRUE, FALSE)           ## whether or not consider Features with complete Channel set
MPpara.arr <-  c("coleff","medres")           ## medianpolish:  MPpara=coleff --> Abundance = overall median + column effect
## medianpolish:  MPpara=medres --> Abundance = overall median + median(residual)
medpolON.arr <-  c("allMix", "indMix")        ## medianpolish to be applied to whole protein block (allMix), or to each Mixture seperately (indMix)
whichSPM.arr <-c("med", "max")                ## for the same PSM per run, mixture and channel
TopNperc <-  5                            ## Top N percent proteins from limma, as insert to heatmap

onlySharedPRT.arr <- c(TRUE, FALSE)           ## Reduce dataset to only shared Proteins between all Mixtures
minPptMIX.arr <- c(1,2)                       ## minimum number of Peptide in each Mixture
mix3rev.arr <- c(TRUE, FALSE)                    


## ensemble mouse
# ens <- useMart("ensembl", "mmusculus_gene_ensembl")

meas <- "new"
PSM2PPT <- FALSE
MPpara <-  "coleff"
medpolON <-  "allMix"




toplist <- data.table(NULL)
fitlist <- data.table(NULL)
key.dt <- data.table(NULL)

# for (ix in 1:length(meas.arr)) {
for (jx in 1:length(FractComb.arr)) {
  for (kx in 1:length(remMix2.arr)) {
    for (lx in 1:length(allChannelsIn.arr)) {
      # for (mx in 1:length(medpolON.arr)) {
      for (nx in 1:length(whichSPM.arr)) {
        for (px in 1:length(onlySharedPRT.arr)) {
          for (qx in 1:length(minPptMIX.arr)) {
            # for(rx in 1:length(mix3rev.arr)) {
              
            # # meas <- meas[ix]
            # FractComb <- FractComb.arr[jx]
            # remMix2 <- remMix2.arr[kx]
            # allChannelsIn <- allChannelsIn.arr[lx]
            # # medpolON <- medpolON[mx]
            # whichSPM <- whichSPM.arr[nx]
            # onlySharedPRT <- onlySharedPRT.arr[px]
            # minPptMIX <- minPptMIX.arr[qx]
            
            print(c(jx, kx, lx, nx, px, qx))

            key.tmp <- data.table(keytag =  paste(c(jx, kx, lx, nx, px, qx), collapse = ".") ,
                                  FractComb = FractComb.arr[jx], 
                                  remMix2 = remMix2.arr[kx],
                                  allChannelsIn = allChannelsIn.arr[lx],
                                  whichSPM = whichSPM.arr[nx],
                                  onlySharedPRT = onlySharedPRT.arr[px],
                                  minPptMIX = minPptMIX.arr[qx],
                                  mix3rev = mix3rev.arr[rx],
                                  meas = meas,
                                  PSM2PPT = PSM2PPT,
                                  MPpara = MPpara,
                                  medpolON = medpolON)
            
            
            key.dt <- rbind(key.dt, key.tmp)
            # "remMix2", "allChannelsIn", "whichSPM", "onlySharedPRT", "minPptMIX")
            
            
            
            ## save toptable and fitlists
            savePhrase <- paste0( "__FracComb_",key.tmp$FractComb,
                                  "__rmNAchnnl_",key.tmp$allChannelsIn,"__whichSPM_",key.tmp$whichSPM,
                                  "__MedPol_",key.tmp$medpolON,"__measurements_",key.tmp$meas,
                                  "__rmMix2_",key.tmp$remMix2, "__onlySharedPRT_",key.tmp$onlySharedPRT,
                                  "__minPptMIX_",key.tmp$minPptMIX)
            
            
            tophits <- fread(file= paste0(getwd(),"/loop_all/tophits_data/topTable_limma", savePhrase, ".csv"),
                              header=TRUE, dec=".", key = "Protein")
            setorder(tophits, Comparison, id)
            tophits$keytag <-  paste(c(jx, kx, lx, nx, px, qx), collapse = ".")
            toplist <- rbind(toplist, tophits)



            limmafit <- fread(file= paste0(getwd(),"/loop_all/statistics_data/fitList_limma", savePhrase, ".csv"),
                             header=TRUE, dec=".", key = "Protein")
            setorder(limmafit, Comparison, id)
            limmafit$keytag <- paste(c(jx, kx, lx, nx, px, qx), collapse = ".")
            fitlist <- rbind(fitlist, limmafit)
            
            
            # savefile <- paste0(getwd(),"/loop_all/volcano_notext/PairwiseVolcanoPlots_",savePhrase,".pdf")
            # volc_cblb(fitList, savefile=savefile, label=FALSE)  
            
          }
        }
      }
    }
  }
}


goi <- data.table(
    Gene= c("Irf4","Nck1","Rela","Cblb") , 
    Protein = c("Q64287", "Q99M51", "Q04207", "Q3TTA7"))


gen.path <- function(tag) {
  key.tmp <- key.dt[keytag %chin% tag,]
  savePhrase <- paste0( "__FracComb_",key.tmp$FractComb,
                        "__rmNAchnnl_",key.tmp$allChannelsIn,"__whichSPM_",key.tmp$whichPSM,
                        "__MedPol_",key.tmp$medpolON,"__measurements_",key.tmp$meas,
                        "__rmMix2_",key.tmp$remMix2, "__onlySharedPRT_",key.tmp$onlySharedPRT,
                        "__minPptMIX_",key.tmp$minPptMIX)
  return(savePhrase)
}

view.volc <- function (tag) {
  savePhrase <- gen.path(tag)
  
  volcano.file <- paste0(getwd(),"/loop_all/volcano/PairwiseVolcanoPlots_",savePhrase,".pdf")
  system(paste0('open "', volcano.file, '"'))
  
  # volcano.notext.file <- paste0(getwd(),"/loop_all/volcano_notext/PairwiseVolcanoPlots_",savePhrase,".pdf")
  # system(paste0('open "', volcano.notext.file, '"'))
}


view.htmp <- function (tag) {
  savePhrase <- gen.path(tag)
  heatmap.file <- paste0(getwd(),"/loop_all/heatmap/TopProt_Heatmap_",savePhrase,".pdf")
  system(paste0('open "', heatmap.file, '"'))
}






tpdt <- toplist[Gene %in% "Cblb" & Comparison == "WT-CBLB"]
tpdt <- tpdt[key.dt, on="keytag"]
ggplot(tpdt, aes(x=allChannelsIn, y=logFC)) + geom_boxplot()



toplist <- toplist[key.dt, on="keytag"]
fitlist <- fitlist[key.dt, on="keytag"]


ggplot(toplist[Comparison == "WT-CBLB" ], aes(x=P.Value, color=allChannelsIn)) +
  geom_histogram(aes(fill=allChannelsIn), alpha=0.5, position="identity", alpha=0.5)


ggplot(fitlist[Comparison == "WT-CBLB" & Gene == "Cblb"], aes(x=logFC, color=allChannelsIn)) +
  geom_histogram(aes(fill=allChannelsIn), alpha=0.5, position="identity", alpha=0.5) 
ggplot(fitlist[Comparison == "WT-CBLB" & Gene == "Cblb"], aes(x=P.Value, color=allChannelsIn)) +
  geom_histogram(aes(fill=allChannelsIn), alpha=0.5, position="identity", alpha=0.5) 





## plot raw protein
work1 <- fread(file = paste0(getwd(),"/loop_all/nonNormdata_work1.csv"),header=TRUE, dec=".", key = "Protein")
dt.raw <- work1[Protein == "Q3TTA7"]
ppt.plot(dt.raw)

## plot normalized protein
work <- fread(file = paste0(getwd(),"/loop_all/Normdata_work.csv"),header=TRUE, dec=".", key = "Protein")
dt.nrm <- work[Protein == "Q3TTA7"]
ppt.plot(dt.nrm)

## plot nonNormalized MedianPolished protein
workf.noNorm <- fread(file = paste0(getwd(),"/loop_all/work_final_noNorm.csv"),header=TRUE, dec=".", key = "Protein")
workf.Norm <- fread(file = paste0(getwd(),"/loop_all/work_final_Norm.csv"),header=TRUE, dec=".", key = "Protein")
dt.raw.mp <- workf.noNorm[Protein=="Q3TTA7"]
dt.nrm.mp <- workf.noNorm[Protein=="Q3TTA7"]
plot.profile(dt.raw, dt.raw.mp)
plot.profile(dt.nrm, dt.nrm.mp)








