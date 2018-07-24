


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
    
    
    
    ####### STEP 2: Peptides, that are used in more than one proteins (Similar to previous step. This is a double check!)
    notUnqPpt <- work0[, c(.(cnt = uniqueN(Protein))), by = list(Peptide)] [cnt!=1, Peptide]
    if (length(notUnqPpt) > 0) work0 <- work0[!Peptide %chin% notUnqPpt,]

    work1 <- copy(work0)
    work1[Abundance=="NaN", Abundance := NA]
    work1 <- char2fact(work1)  
    
    
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
      # table(cnt.wch0$cnt)
      work2 <- rbindlist(list(wch0[,-c("cnt")], wch2[,-c("cnt", "mean_int")] ))
      
    } else { work2 <- copy(work1) }
    
    
    ########################################################
    ################################### remove Features with only NAs in all channels
    cnt.work2 <- work2[, c(.(cnt = sum(!is.na(Abundance)))), by = list(Mixture, Protein, Feature)] 
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
    work <- VSNnorm(work3)
        setnames(work, old="Abundance", new="Abundance.noNorm")
        setnames(work, old="Abundance.norm", new="Abundance")
        work <- char2fact(work)
    ################################################################# 
    
    
    
    ################################################################# 
    ## Median Polish 
    # median polish over all mixture at once
    if (medpolON == "allMix") work.mp <- MedPolAll(work)
    # median polish for each mixture separetly
    if (medpolON == "indMix") work.mp <- MedPolInd(work) 
    ################################################################# 
    
    
    
    #################################################################   
    ## final work  
    workf <- work.mp[
      unique(work[, list(Protein, Mixture, Channel, Condition, BioReplicate)]),
      on=c("Protein", "Mixture", "Channel")]
    
    # tmp <- data.table(Peptide = as.character(unique(work$Peptide)), PeptideN=NA)
    # tmp[, PeptideN := sapply(Peptide, function(x) .spltpep(x,n=9))]
    # work <- work[tmp, on="Peptide"]
    
    
    ################################################################# 
    ## extract bioloical parameters  ..> genes
    # ens <- useMart("ensembl", "hsapiens_gene_ensembl")
    trans = getBM(attributes = c("uniprotswissprot", "mgi_symbol", "ensembl_gene_id", #"ucsc", "mgi_symbol"
                                 "external_gene_name",  "entrezgene", "description"), 
                  # filters = "uniprotswissprot", values = "Q64213", 
                  filters = "uniprotswissprot", values = unique(workf$Protein),
                  uniqueRows=FALSE, mart = ens)
    trans <- as.data.table(trans)
    
    
    setkey(workf, "Protein")
    setkey(trans, "uniprotswissprot")
    workf[trans, Gene := external_gene_name]
    # workf[Protein %in% "Q8CDM1" , Gene := "Atad2"]
    ################################################################# 
    
    ## test median polish -- plots
    # dd <- work[Protein == "Q6ZQK0"]
    # ddd <- workf[Protein == "Q6ZQK0"]
    # plot.profile(dd,ddd)
    # dd <- work[Protein == "Q9QZ67"]
    # ddd <- workf[Protein == "Q9QZ67"]
    # plot.profile(dd,ddd)             
    
    #
    # dd <- work[Protein=="Q9D8V7"]
    # ddd <- workf[Protein=="Q9D8V7"]
    # 
    # dd <- work[Protein=="Q6P2L6"]
    # ddd <- workf[Protein=="Q6P2L6"]
    # 
    # dd <- work[Protein=="Q64287"]
    # ddd <- workf[Protein=="Q64287"]    
    # plot.profile(dd,ddd)
    # 
    # ddt <- work[Protein=="Q64287"]
    # prt.plot(ddt)
    
    # dd <- work[Protein == "Q3TTA7"]
    # ddd <- workf[Protein == "Q3TTA7"]
    # plot.profile(dd,ddd)
    
    
    
    
    ################################################################# 
    ## replace wt and cblb tags in mix2
    if (mix3rev) {
      
      workf[Mixture ==2 & Condition == "WT.0h", Condition := "__Cblb.0h"]
      workf[Mixture ==2 & Condition == "WT.24h", Condition := "__Cblb.24h"]
      workf[Mixture ==2 & Condition == "WT.48h", Condition := "__Cblb.48h"]
      
      workf[Mixture ==2 & Condition == "Cblb.0h", Condition := "WT.0h"]
      workf[Mixture ==2 & Condition == "Cblb.24h", Condition := "WT.24h"]
      workf[Mixture ==2 & Condition == "Cblb.48h", Condition := "WT.48h"]
      
      workf$Condition <- gsub("__", "" ,workf$Condition) 
      workf$Condition <- as.factor(workf$Condition)
      
    }
    ################################################################# 
    
    
    
    
    ## this will be needed  
    workf0 <- copy(workf)
    ################################################################# 
    ## statistical analysis -- limma
    ## TopNperc = top N percent hits in topList ... default TopNperc = 5
    limmaOut <- statSig(workf, TopNperc)
    topList <- limmaOut[[1]]
    fitList <- limmaOut[[2]]
    ################################################################# 
    
    
    
    
    ################################################################# 
    ## Generate PCA plots
    prtvar <- workf0[, c(.(var = var(Abundance))), by = list(Protein)]
      tmp <- prtvar[order(-var)]
      tmp <- tmp[1:500,]
      tmp <- workf0[Protein %in% tmp$Protein,]
    
        wk <- data.table::dcast(tmp, Protein ~ Condition + Mixture + BioReplicate, value.var = "Abundance")
        wk$na_count <- apply(wk[,-1], 1, function(x) sum(is.na(x)))
        wk <- wk[na_count ==0]
        mat <- as.matrix(wk[,-c("Protein","na_count")])
        pcaF = prcomp(t(mat), scale. = FALSE, center = TRUE)
    
    ## plot pca     
    pcaplot <- pca.plot(pcaF)
    ################################################################# 
    
    
    
    
    ################################################################# 
    ## Generate Heatmaps            
    ## get the top N Proteins of each comparison
    topN <- 10
    topPrtCom = topList[topList[, .I[1:topN], Comparison]$V1, list(Comparison, Protein, Gene)]   
    hm.protein <- unique(topPrtCom$Protein)
    
    reftb <- unique(workf0[,c("BioReplicate", "Condition", "Mixture", "Channel")])
    
      hm.dt <- workf[Protein %chin% hm.protein] 
      if (!"Cblb" %in% hm.dt$Gene) { 
        hm.dt <- rbind(workf[Gene %chin% c("Irf4","Cblb")], hm.dt) ## add cblb and Irf4
      }
      
      hm.dt[, Rep := tstrsplit(BioReplicate, "_",  fixed=TRUE)[3]]
      hm.dt[, xTag := do.call(paste, c(.SD, sep = "\nRep")), .SDcols=c("Condition", "Rep")]
      hm.dt[, yTag := do.call(paste, c(.SD, sep = " _ ")), .SDcols=c("Protein", "Gene")]
      
      hm.dt.w <- data.table::dcast(hm.dt, yTag ~ xTag, value.var="Abundance" )
      hm.mat <- as.matrix(hm.dt.w[,-1])
      rownames(hm.mat) <- hm.dt.w$yTag
    
    
    ## plot heatmaps
    hmOut <- hmap.plot(hm.mat)
    
    ################################################################# 
    
    
    
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
    
    
    # pdf(paste0(getwd(),"/loop_all/heatmap/TopProt_Heatmap_",savePhrase,".pdf"), width = 12, height = 11)
    # print(hm.plot)
    # dev.off()
    # 
    # pdf(paste0(getwd(),"/loop_all/heatmap_allclust/TopProt_Heatmap_",savePhrase,".pdf"), width = 12, height = 11)
    # print(hmallclust.plot)
    # dev.off()
    
    pdf(paste0(getwd(),"/loop_all/heatmap_all/TopProt_Heatmap_",savePhrase,".pdf"), width = 17, height = 13)
    draw(hmOut, annotation_legend_side = "bottom")
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

