


# rm(list=ls())
setwd("/Users/farhad/_Rspace/_prj_Wolf_Isabelle/")
ddir <- "/Users/farhad/_Rspace/_prj_Wolf_Isabelle/data/"
adir <- "/Users/farhad/_Rspace/_prj_Wolf_Isabelle/analysis/"
rdir <- "/Users/farhad/_Rspace/_prj_Wolf_Isabelle/scripts/"

source(file="scripts/libs.R")
source(file="scripts/tmtPlots.R")
source(file="scripts/tmtFuncs.R")
source(file="scripts/statSig.R")

source(file="scripts/tmtReform.R")




## Dynamic parameters for version chack
###################################
# meas.arr <- c("old" ,"new")                   ## Which measurement
betweenFracNorm.arr <- c("median", "quantile", "none")
betweenChnNorm.arr = c("median","vsn", "qntVSN.noCal", "medVSN.noCal")


Frac10.arr = c(TRUE, FALSE) 
PSM2PPT.arr = c(FALSE)                  ## PSM to Peptide (Ratio > Median) 
FractComb.arr <-  c("max", "sum", "single")      ## single with VSN::Strata and VSN, max after MedianEq per Mix + Fract
remMix2.arr <- c(FALSE, TRUE)                 ## remove mixture 2 or not (high missingness in mix2)
allChannelsIn.arr <- c(TRUE, FALSE)           ## whether or not consider Features with complete Channel set
MPpara.arr <-  c("coleff")           ## medianpolish:  MPpara=coleff --> Abundance = overall median + column effect
## medianpolish:  MPpara=medres --> Abundance = overall median + median(residual)
medpolON.arr <- c("allMix")        ## MP to be applied to all mixtures or to each Mixture seperately 
whichPSM.arr <- c("med")                ## for the same PSM per run, mixture and channel
TopNperc <-  5                                ## Limma TopTable: top N percent hits

onlySharedPRT.arr <- c(TRUE, FALSE)           ## Reduce dataset to only shared Proteins between all Mixtures
onlySharedFEAT <- c(TRUE, FALSE)

minPptMIX.arr <- c(1,2)                       ## minimum number of Peptide in each Mixture
# mix3rev.arr <- c(TRUE)                    




# Constatnt parameters for version chack
betweenFracNorm = "none"
betweenChnNorm = "vsn"
vsn.Calib <- "affine"
Frac10 = c(TRUE, FALSE)
PSM2PPT <- FALSE
MPpara <-  "coleff"
medpolON <-  "allMix"
allChannelsIn <- FALSE
whichPSM <- "med"
# remMix2 <- FALSE      
minPptMIX <- 0
# FractComb <- "max"
onlySharedPRT <- FALSE
onlySharedFEAT <- FALSE
TopNperc <-  5 


key.dt <- data.table(NULL)

for (ix in 1:length(Frac10.arr)) {
for (jx in 1:length(FractComb.arr)) {
for (kx in 1:length(allChannelsIn.arr)) {
for (lx in 1:length(remMix2.arr)) {
# for (mx in 1:length(onlySharedPRT.arr)) {
for (nx in 1:length(minPptMIX.arr)) {
for (px in 1:length(betweenFracNorm.arr)) {
for (qx in 1:length(betweenChnNorm.arr)) {
        
    
    Frac10 <- Frac10.arr[ix]
    FractComb <- FractComb.arr[jx]
    allChannelsIn <- allChannelsIn.arr[kx]
    remMix2 <- remMix2.arr[lx]
    # onlySharedPRT <- onlySharedPRT.arr[mx]
    minPptMIX <- minPptMIX.arr[nx]
    betweenFracNorm <- betweenFracNorm.arr[px]
    betweenChnNorm <- betweenChnNorm.arr[qx]
    
    
    keytag <- paste(c(ix, jx, kx, lx, nx, px, qx), collapse = ".") 
    key.tmp <- data.table(keytag = keytag,
                          Frac10 = Frac10[ix],
                          FractComb = FractComb.arr[jx],
                          allChannelsIn = allChannelsIn.arr[kx],
                          remMix2 = remMix2.arr[lx],
                          # onlySharedPRT = onlySharedPRT.arr[mx],
                          minPptMIX = minPptMIX.arr[nx],
                          betweenFracNorm = betweenFracNorm.arr[px],
                          betweenChnNorm = betweenChnNorm.arr[qx]
                          )
    
 
    savePhrase <- paste0( keytag,"__Frac10_",Frac10,"__Frac_",FractComb,
                          "__rmNAch_",allChannelsIn,
                          "__rmM2_",remMix2,
                          "__minPpt_",minPptMIX, "__FracNrm_", betweenFracNorm, 
                          "__ChnNrm_", betweenChnNorm) 
    
    savePhrase <- gsub("TRUE","T", savePhrase)
    savePhrase <- gsub("FALSE","F", savePhrase)
    key.dt <- rbind(key.dt, key.tmp)
    
    print(key.tmp)
  
    


    ## initiating work data set from PSM data set
    work0 <- psm.ready[, list(Frac.i, Run, Channel, Protein, Peptide, Feature, 
                          Charge, Proptide, Condition, BioReplicate, 
                          Mixture, Abundance0, Gene)]
    
    
    ########################################################
    ## remove fractions 1 & 10 from all mixtures       
    if (Frac10) {
    RunOut <- c("1.2", "1.3", "1.4", "10.2", "10.3", "10.4") 
    work0 <- work0[!Run %in% RunOut]
    }
    
    
    ########################################################
    ## Filterring procedure
    
    ####### STEP 2: Peptides, that are used in more than one proteins (Similar to previous step. This is a double check!)
    notUnqPpt <- work0[, c(.(cnt = uniqueN(Protein))), by = list(Peptide)] [cnt!=1, Peptide]
    if (length(notUnqPpt) > 0) work0 <- work0[!Peptide %chin% notUnqPpt,]

    work1 <- copy(work0)
    work1[Abundance0=="NaN", Abundance0 := NA]
    work1 <- char2fact(work1)  
    
    
    ## NA count per feature and Run -- before Fraction combintaion
    NAcnt.work1 <- work1[, c(.(NAper = sapply( sum(is.na(Abundance0)) / length(Abundance0) , function(x) round(x,3) ) ) ,
                             .(NAcnt = sum(is.na(Abundance0))),
                             .(TTLcnt = length(Abundance0))) ,
                         by = list(Feature)]
    # table(NAcnt.work1$NAcnt)
    # table(NAcnt.work1$TTLcnt)
        NA.tb <- as.data.table(table(NAcnt.work1$NAper))
        NAfeat <- round(100 * (sum(NA.tb[V1!=0, N]) / sum(NA.tb$N)) ,3) # percentage of Features with at least 1 NA
        mnNA <- round(100 * (summary(as.numeric(NA.tb[V1!=0, V1]))["Min."]),4)
        mxNA <- round(100 * (summary(as.numeric(NA.tb[V1!=0, V1]))["Max."]),4)
        
        print(paste0(NAfeat, "% of the Features before fraction combination have at least 1 NA value across all Runs and Channels."))
        print(paste0("Degree of missingness in Features varies from ~",mnNA, "% to ~", mxNA, "%"))
    
    
    
    ########################################################
    ## Dealing with multiple PSM measurements within and across fractions.
    ## Below: Number of measurements (intensity values) per Feature per Run. This must ideally be equal to number of Channels (here 6)
    ## We have to decide what to do when a Feature is measured more than once. 
    
    
    ###############################
    ##STEP 1: mult. measurements for same Feature, same Mixture, same Channel and same Run: calc. Mean OR Median
    if (remMix2) { work1 <-  work1[Mixture != 2] } # remove mixture 2
    # work11 <- copy(work1)
    
    if (toupper(whichPSM) == "MEAN") {
      work1[, Abundance0 := mean(Abundance0, na.rm=TRUE), by=list(Feature, Mixture, Run, Channel)] #Mean
    }
    if (toupper(whichPSM) == "MED") {
      work1[ ,Abundance0 := median(Abundance0, na.rm=TRUE), by=list(Feature, Mixture, Run, Channel)] #Median
    }
    if (toupper(whichPSM) == "MAX") {        
      work1[ ,Abundance0 := max(Abundance0, na.rm=TRUE), by=list(Feature, Mixture, Run, Channel)] #max
    } 
    work1 <- unique(work1)
    

    ###############################
    ####### STEP 2: From PSM to Peptide level:
    ###############################
    ## Number of Charge states per Peptide per Protein
    cnt.work1 <- work1[, .(n = uniqueN(Charge)), by = list(Proptide)] 
    
    if (PSM2PPT) {
      ## CASE 1: Peptides with 2 PSM > i. compute PSM ratios (in log scale: ratio = PSM1 - PSM2)
      ##                               ii. compute median of ratios (medRatio)
      ##                               iii. when PSM2 is NA & PSM1 nonNA: subtract medRatio from PSM1 (PSM2 = PSM1 - medRatio)
      ##                               iv . when PSM1 is NA & PSM2 nonNA: add medRatio to PSM2 (PSM1 = PSM2 + medRatio)
      
      ## subset data to Peptides with 2 PSMs
      pptpsm2 <- cnt.work1[n == 2, Proptide] # Peptides with 2 PSMs
      work.2psm <- work1[Proptide %in% pptpsm2, ] # corresponding data block
      work.2psm.w <- data.table::dcast(work.2psm, Proptide+Feature ~ Mixture + Run + Channel, value.var = "Abundance0") #to wide format
      work.2psm.merged <- work.2psm.w[, f.psm2ppt(.SD) , by = .(Proptide)  ]  
      
      ## merge with original data and output
      work1 <- merge(work.2psm.merged, work.2psm[,c("Feature", "Mixture", "Run", "Channel", "Abundance0")],
                     by = c("Feature", "Mixture", "Run", "Channel"), all.x = TRUE)
      
    } # end of PSM2PPT
    
    
    
    ########################################################
    ################################### Median & Quantile normalization between fractions and mixtures
    work1. <- copy(work1)
    if (betweenFracNorm != "none") {
      
      if (toupper(betweenFracNorm) == "MEDIAN") {
          ## median Norm. across Fractions 
          medFrac <- work1.[, .(medFrac = median(Abundance0, na.rm=TRUE)), by = c("Run", "Channel")]
          work1. <- merge(work1., medFrac, by= c("Run", "Channel"), all.x=TRUE)
          work1.[, Abundance := Abundance0 - medFrac + median(medFrac, na.rm=TRUE)]
          
          FractComb <- "max" ## by median normalization, fraction combination must be Max 
      }

      if (toupper(betweenFracNorm) == "QUANTILE") {
        ## quantile Norm. across Fractions & Channels 
          ## !!!!! Check the quantile normalization for NA hendling. 02.08.2018
          tmp.dt <- data.table::dcast(work1., Protein + Feature ~ Channel + Run, value.var = "Abundance0")
          tmp.mat <- as.matrix(tmp.dt[,!1:2])
          
            norm.mat = normalize.quantiles(tmp.mat)
            colnames(norm.mat) <- colnames(tmp.mat)
              norm.dt.w <- cbind(tmp.dt[,1:2], as.data.table(norm.mat))
              norm.dt <- melt.data.table(norm.dt.w, id.vars = c("Protein","Feature"), 
                                         variable.name = "Run", value.name = "Abundance")
              norm.dt[, c("Channel","Run") := tstrsplit(Run, "_",  fixed=TRUE)]
              norm.dt <- char2fact(norm.dt)
          work1. <- merge(work1., norm.dt, by= c("Protein","Feature", "Run", "Channel"), all.x=TRUE)              
      }      
              
      
    } else {
      work1.[, Abundance := Abundance0]
    }
    
    work1. <- work1.[,-c("medFrac")]
    
    
    ########################################################
    ################################### work summary before combination    
    ## NAs per condition, mixture and rep
    # NAcond1 <- work1.[, c( .(NA.N = sum(is.na(Abundance)) ),
    #                       .(NA.Perc = length(Abundance)) ), 
    #                  by = list(Mixture, BioReplicate, Condition)] 
    # NAcond1[, NA.Perc := round(NA.N/NA.Perc, 4)*100 ]
    # NAcond1[, FracComb. := "Frac. not Combined" ]
    
    
    
    ########################################################
    ################################### Fraction combination
      work2 <- copy(work1.)
      sharedPSM.Mix <- work2[, c(.(cnt = uniqueN(Run))), by = list(Feature, Mixture)] #table(sharedPSM.Mix$cnt)
            
            # ###### plots
            # ## distribution of PSMs over fractions -- plot
            # sharedPSM.dt <- as.data.table(table(sharedPSM.Mix$cnt))
            # sharedPSM.dt[, perc := round(N/sum(N),4)*100 ]
            # sharedPSM.dt[, txt := paste0("#",N,"\n" ,perc,"%")]
            # 
            # MultFrac <- 
            #   ggplot(sharedPSM.dt, aes(x=as.integer(V1), y=N)) + 
            #   geom_bar(stat= "identity") + 
            #   geom_text(col="red", vjust=-0.3, aes(label=txt)) +
            #   scale_x_discrete(name ="Fractions count", limits=sharedPSM.dt$V1) +
            #   labs(title = "How many PSMs measured in how many Fractions!", y = "PSM count")  +
            #   scale_y_continuous(limits=c(0,85000))
            # ######
            
            
    if (toupper(FractComb) != "NONE") {
      # check: ppt.plot(work[Proptide == "P52480_lAPITSDPTEAAAVGAVEASFk"])
      wch0 <- work2[sharedPSM.Mix[cnt==1,], on=c("Feature", "Mixture")][,-"cnt"] #unique measurements
            
      if (toupper(FractComb) == "SINGLE") { # continue with Proteins with only 1 Fractions
        work2 <- copy(wch0)
      } else {
        
        wch  <- work2[sharedPSM.Mix[cnt>1,], on=c("Feature", "Mixture")][,-"cnt"]  #redundant measurements
        
            if (toupper(FractComb) != "SUM") {
              
                wch1 <- wch[ , .(mean_int = mean(Abundance, na.rm=TRUE)), by = c("Protein","Feature","Mixture","Run")] 
                  
                  if (toupper(FractComb) == "MAX") { #get the maximum fraction
                    wch.frc  <- wch1[wch1[, .I[which.max(mean_int)],  by=c("Protein","Feature","Mixture")]$V1]
                  }
                  
                  if (toupper(FractComb) == "MED") { #get the median fraction ... see: which.medX
                    wch.frc <- wch1[wch1[, .I[which.medX(mean_int)], by=c("Protein","Feature","Mixture")]$V1]
                  }
                
                wch2 <- wch[wch.frc, on=c("Protein","Feature","Mixture","Run")] 
                cnt.wch2 <- wch2[, c(.(cnt = sum(!is.na(Abundance)))), by = list(Mixture, Protein, Feature)] 
                cnt.wch0 <- wch0[, c(.(cnt = sum(!is.na(Abundance)))), by = list(Mixture, Protein, Feature)] 
                
                work2 <- rbindlist(list(wch0, wch2[,-c("cnt", "mean_int")] ))
              
            
                ######
                ## distribution of PSMs over fractions -- plot
                MaxFrac.dt <- as.data.table(table(work2$Run)/6)
                MaxFrac.dt[, c("Run","Mixture") := tstrsplit(V1, ".",  fixed=TRUE)]
                MaxFrac.dt <- MaxFrac.dt[order(as.integer(V1), Mixture)]
                
                MaxFrac <- 
                  ggplot(MaxFrac.dt, aes(x=Run, y=N, label=N)) + 
                  facet_grid(~Mixture) +
                  geom_bar(stat= "identity") + 
                  geom_text(col="red", vjust=0.5, hjust=-0.1, angle = 90) +
                  scale_x_discrete(name ="Run ID", limits=as.factor(MaxFrac.dt$Run)) +
                  labs(title = "How many Max PSMs in each Run", y = "maxPSM count")  
                ######
                
            } ## end non Sum
        
        if (toupper(FractComb) == "SUM") { 
          
          wch.sum <- unique(wch[, c(.(Abundance = sum(Abundance))), by = c("Channel", "Mixture", "Feature", "Protein")])
          wch <- unique(wch[, -c("Run", "Frac.i", "Abundance", "Abundance0")])
          wch0 <- wch0[, -c("Run", "Frac.i")]
          wch.sum <- unique(wch[wch.sum, on=c("Protein","Feature","Mixture", "Channel")] )
          
          work2 <- rbind(wch0[,-c("Abundance0")], wch.sum )
        }  ## end sum
        
      } # end FractComb: SINGLE
    
    } else { 
      work2 <- copy(work1) 
    }
    
    work2 <- work2[,-c("Abundance0", "Run", "Frac.i")]
    
    ########################################################
    ################################### remove Features with only NAs in all channels
    cnt.work2 <- work2[, c(.(cnt = sum(!is.na(Abundance)))), by = list(Mixture, Protein, Feature)] 
    work3 <- work2[cnt.work2[cnt != 0,-"cnt"], on=c("Protein", "Mixture", "Feature")]
    
    
    
            ########################################################
            ################################### work summary after combination
            ## NAs per condition, mixture and rep
            # NAcond3 <- work3[ , c( .(NA.N = sum(is.na(Abundance))),
            #                        .(NA.Perc = round(sum(is.na(Abundance))/uniqueN(Feature), 4)*100 ),
            #                        .(cnt.Peptide = uniqueN(Peptide)),
            #                        .(cnt.Protein = uniqueN(Protein)),
            #                        .(cnt.Feature = uniqueN(Feature)) ),
            #                   by = list(Mixture, BioReplicate, Condition)] 
            # NAcond3[, FracComb. := "Frac. Combined" ]
            # NAcond.dt <- data.table::melt(NAcond3, id.vars =  c("Mixture", "BioReplicate", "Condition","FracComb."),
            #                               variable.name="Para", 
            #                               value.name="Val")
            # 
            # 
            # perCondPlot <-
            #   ggplot(NAcond.dt, aes(x=Condition, y=Val, fill=Mixture, label=Val)) + 
            #   facet_grid(Para ~ BioReplicate, scales = "free") +
            #   geom_bar(stat= "identity") + 
            #   geom_text(col="black", vjust=0.6, hjust=0.8, angle = 45, size=3.5) +
            #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
            #   scale_fill_brewer(palette = "Paired") 
            ######
    
    
    
    ########################################################
    ################################### remove Features with not complete channel-set

      cnt.work3 <- work3[, c(.(cnt = sum(!is.na(Abundance)))), by = list(Mixture, Protein, Feature)] 
      
                # fullCHN <- as.data.table(table(cnt.work3[, list(Mixture,cnt)])) ## missingness per feature per mixture
                # fullCHN[, ttl := sum(N), by=Mixture ]
                # fullCHN[, Percentage := round(100*N/ttl,2) ]
                # fullCHN[, txt := paste0("#",N,"\n" ,perc,"%")]
                # 
                # MultChnnl <- 
                #   ggplot(fullCHN, aes(x=cnt, y=Percentage, fill = Mixture)) + 
                #   geom_bar(stat= "identity", position="dodge") + 
                #   geom_text(vjust=0.3, aes(label=Percentage), size=4, angle=90, position = position_dodge(width = 1)) +
                #   scale_x_discrete(name ="Channel-count", limits=fullCHN$V1) +
                #   scale_fill_brewer(palette = "Paired") +
                #   labs(title = paste0(
                #     "Count/Percentage of Features vs. number of non NA measurements per Mixture.\n",
                #     "Value of 6 means that all Channels of a given Feature have non NA values."), 
                #     y = "PSM count")  
      
    if (allChannelsIn) {
      work3 <- work3[cnt.work3[cnt == 6, -"cnt"], on=c("Protein", "Mixture", "Feature")]
    }
    
    
    ########################################################
    ################################### continue with Proteins exist in all Mixrture    
    if (onlySharedPRT) {
      nPrtMix <- work3[, c(.(cnt = uniqueN(Mixture))), by = list(Protein)]
      work3 <- work3[Protein %in% nPrtMix[cnt==uniqueN(work3$Mixture), Protein] ]
      # work3 <- work3[nPrtMix[cnt==uniqueN(work3$Mixture), -"cnt"], on="Protein"]
    }
    

      
    ########################################################
    ################################### continue with Proteins exist in all Mixrture    
    if (onlySharedFEAT) {
      nFeatMix <- work3[, c(.(cnt = uniqueN(Mixture))), by = list(Feature)]
      work3 <- work3[Feature %in% nFeatMix[cnt==uniqueN(work3$Mixture), Feature] ]
      # work3 <- work3[nPrtMix[cnt==uniqueN(work3$Mixture), -"cnt"], on="Protein"]
    }
      
      
      
    ########################################################
    ################################### continue with Proteins with a mimimum number of Peptides in Mixture (minPptMIX)   
    if (minPptMIX != 0) {
      nPptMix <- work3[, c(.(cnt = uniqueN(Peptide))), by = list(Protein, Mixture)]
      nPptMix.w <- dcast(nPptMix, Protein ~ Mixture, value.var = "cnt")
      nPptMix.w[is.na(nPptMix.w)] <- 0
         nPptMix.l <- data.table::melt(nPptMix.w, id.vars =  "Protein",
                                    variable.name="Mixture", 
                                    value.name="cnt")
      
      prtOUT <- unique(nPptMix.l[cnt < minPptMIX, Protein])
      if (length(prtOUT) != 0 ) work3 <- work3[!Protein %in% prtOUT]
    } 
    

      
      # ################################################################# 
      # ## NA count per condition per Mixture   
      # work3.w  <- data.table::dcast(work3, Protein + Peptide + Charge ~ Mixture + Channel, value.var = "Abundance")
      # cols.ch <- c("Protein","Peptide","Charge")
      # cols <- names(work3.w[,-cols.ch, with=FALSE])
      # ## number of NAs per column and row
      # na.col <- work3.w[, lapply(.SD, function(x) sum(is.na(x))), .SDcols = cols]       #NA per col    table(na.col)
      # na.col <- melt(na.col, variable.name = "Run", value.name = "NA.count", verbose = FALSE)
      # na.col[, c("Mixture", "Channel") := tstrsplit(Run, "_",  fixed=TRUE)]
      # na.col <- na.col[unique(work3[,c("Mixture","Channel","Condition", "BioReplicate")]) , on=c("Mixture","Channel")]
      # 
      # NAcntBar <-
      #   ggplot(na.col, aes(x=Condition, y=NA.count, col=Condition)) +
      #   scale_color_brewer(palette = "Dark2") +
      #   theme_bw() +
      #   # geom_bar(aes(y = (..count..)/sum(..count..)))
      #   geom_bar(stat="identity",position = "dodge2") +
      #   facet_wrap(~Mixture, scales = "free_x") +
      #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position="right") 
      
      
      
      
    ################################ 
    ################################################################# 
      work3. <- copy(work3)
      
      if ( toupper(betweenChnNorm) == "VSN") { ## VSN normalization  
        
        vsnOut <- VSNnorm(work3., calib = "affine") ## vsn.Calib = "affine" , "none"
            work    <- vsnOut[[1]] ## dataset with vsn abundance
            vsnTest <- vsnOut[[2]] ## rank mean v. SD from VSN
            work <- char2fact(work)
      }  
        
      
      if ( toupper(betweenChnNorm) == "QNTVSN.NOCAL") { ## between channel Quantile + VSN with calib="none"
        ## quantile between channels
        tmp.dt <- data.table::dcast(work3., Protein + Feature  ~ Channel + Mixture, value.var = "Abundance")
        tmp.mat <- as.matrix(tmp.dt[,!1:2])
        norm.mat = normalize.quantiles(tmp.mat)
        colnames(norm.mat) <- colnames(tmp.mat)
        norm.dt.w <- cbind(tmp.dt[,1:2], as.data.table(norm.mat))
        norm.dt <- melt.data.table(norm.dt.w, id.vars = c("Protein","Feature"), 
                                   variable.name = "chn", value.name = "Abundance")
        norm.dt[, c("Channel", "Mixture") := tstrsplit(chn, "_",  fixed=TRUE)]
        norm.dt[,chn := NULL]
        work3. <- merge(work3.[,-c("Abundance")], norm.dt, by= c("Protein","Feature", "Mixture", "Channel"), all.x=TRUE)
        
        ## VSN between channels
        vsnOut <- VSNnorm(work3., calib = "none") ## calib = "none" , already with quantile
        work    <- vsnOut[[1]] ## dataset with vsn abundance
        vsnTest <- vsnOut[[2]] ## rank mean v. SD from VSN
        work <- char2fact(work)
      }
      
      
      if ( toupper(betweenChnNorm) == "MEDVSN.NOCAL") { ## median normalization
        
        medFrac <- work3.[, .(medFrac = median(Abundance, na.rm=TRUE)), by = c("Mixture", "Channel")]
        work3. <- merge(work3., medFrac, by= c("Mixture", "Channel"), all.x=TRUE)
        work3.[, Abundance := Abundance - medFrac + median(medFrac, na.rm=TRUE)]
        
        ## VSN between channels
        vsnOut <- VSNnorm(work3., calib = "none") ## calib = "none" , already with quantile
        work    <- vsnOut[[1]] ## dataset with vsn abundance
        vsnTest <- vsnOut[[2]] ## rank mean v. SD from VSN
        work <- char2fact(work)
      }
      
      
        
      if ( toupper(betweenChnNorm) == "MEDIAN") { ## median normalization
        
        medFrac <- work3.[, .(medFrac = median(Abundance, na.rm=TRUE)), by = c("Mixture", "Channel")]
        work <- merge(work3., medFrac, by= c("Mixture", "Channel"), all.x=TRUE)
        work[, i.Abundance := Abundance]
        work[, Abundance := Abundance - medFrac + median(medFrac, na.rm=TRUE)]
        
        mat.w <- dcast(work, Protein + Peptide + Charge ~ Mixture + Channel, value.var = "Abundance")
        mat <- as.matrix(mat.w[,-c("Protein","Peptide","Charge")]) 
        vsnTest <- meanSdPlot(mat, plot=FALSE)$gg + theme_bw() +  scale_fill_distiller(palette = "Spectral") 
      }     
      
      
    ################################################################# 
       
      sctPlot <- ## scatterpolt median v. quantile
        ggplot(work , aes(Abundance, i.Abundance)) + 
        stat_binhex(aes(alpha=..count.., fill=Condition)) +
        geom_abline(intercept = 0, slope = 1, col="darkred", alpha=0.6) +
        # geom_smooth(se = FALSE, method = lm, alpha=0.2, size=0.5, aes(col=Condition)) +
        theme_bw() + facet_wrap(~Mixture) + coord_fixed()
      
        
      pDen.chnMed <- pDen(work, x="Abundance", col="Condition")
      pBox.chnMed <- pBox(work, x="Condition", y="Abundance", col="BioReplicate")  
      
      
    ################################################################# 
    ## Median Polish 
    # median polish over all mixture at once
    if (medpolON == "allMix") {
      work.mp <- MedPolAll(work, MPpara="coleff")
    }
    # # median polish for each mixture separetly
    # if (medpolON == "indMix") {
    #   work.mp <- MedPolInd(work)[[1]]
    # }
    ################################################################# 
    
    
    #################################################################   
    ## final work  
    workf <- work.mp[
      unique(work[, list(Gene, Protein, Mixture, Channel, Condition, BioReplicate)]),
      on=c("Protein", "Mixture", "Channel")]


          
      # test some proteins
      ################################################################# 
      # dd <- work[Protein=="Q64287"]
      # ddd <- workf[Protein=="Q64287"]
      # plot.profile(dd,ddd)
      # # 
      # dd <- work[Protein == "Q3TTA7"]
      # ddd <- workf[Protein == "Q3TTA7"]
      # plot.profile(dd,ddd)
      
      
      
      
    ## this will be needed  
    workf0 <- copy(workf)
    ################################################################# 
    ## statistical analysis -- limma
    ## TopNperc = top N percent hits in topList ... default TopNperc = 5
    limmaOut <- statSig(workf, TopNperc=5)
    topList <- limmaOut[[1]]
    fitList <- limmaOut[[2]]
    
    ## plot volcano
    volc.all <- volc.p.all(fitList)
    volc.ind <- volc.p.ind(fitList, top=topList)
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
    topN <- 20
    topPrtCom = topList[topList[, .I[1:topN], Comparison]$V1, list(Comparison, Protein, Gene)]  
    
    
        # #####################################################
        # ## Generating two fixed toplists of proteins for heatmaps.
        # ## The two lists contain top 100 proteins from two comparisons: Cblb.24h-WT.24h & WT-CBLB
        # ## The version used for generating toplists has following parameters:
        # ## "__mix3rev___FracComb_max__rmNAchnnl_FALSE__whichPSM_med ...
        # ## ...   __MedPol_allMix__measurements_new__rmMix2_FALSE__onlySharedPRT_FALSE__minPptMIX_1"
        #      prtFIX.WtKo24 <- data.table(Protein=topList[Comparison == "Cblb.24h-WT.24h", Protein][1:100])
        #      prtFIX.WtKo   <- data.table(Protein=topList[Comparison == "WT-CBLB", Protein][1:100])
        #       write.csv(prtFIX.WtKo24, file = paste0(adir,"/loop_all/prtFIX_WtKo24_FracMax_minPep1.csv"))
        #       write.csv(prtFIX.WtKo,   file = paste0(adir,"/loop_all/prtFIX_WtKo_FracMax_minPep1.csv"))
        # #####################################################
    
    ## Reading Proteins Fixed toplists  ... fsetdiff(prtFIX.WtKo24,prtFIX.WtKo) ... 60 different Proteins
    ## fixed Protein list from wt24-cblb24
    prtFIX.WtKo24 <- fread(file= paste0(adir,"/loop_all/prtFIX_WtKo24_FracMax_minPep1.csv"), header=TRUE)$Protein
    ## fixed Protein list from wt-cblb
    prtFIX.WtKo <- fread(file= paste0(adir,"/loop_all/prtFIX_WtKo_FracMax_minPep1.csv"), header=TRUE)$Protein
    
    ##Protein list from all comparisons
    hm.proteinAll <- unique(topPrtCom$Protein)
    ##Protein list from all wt v. cblb comparisons
    hm.proteinWtKo <- unique(topPrtCom[Comparison %like% "Cblb" & Comparison %like% "WT", Protein])
    
    
    
    heatmap.mat <- function(prtl) {
    reftb <- unique(workf0[,c("BioReplicate", "Condition", "Mixture", "Channel")])
    
      hm.dt <- workf[Protein %chin% prtl] 
      if (!"Cblb" %in% hm.dt$Gene) { 
        hm.dt <- rbind(workf[Gene %chin% c("Cblb")], hm.dt) ## add cblb and Irf4
        hm.dt[Gene == "Cblb", Gene := "**Cblb**"]
      }
    
      if (!"Irf4" %in% hm.dt$Gene) { 
        hm.dt <- rbind(workf[Gene %chin% c("Irf4")], hm.dt) ## add cblb and Irf4
        hm.dt[Gene == "Irf4", Gene := "**Irf4**"]
      }
    
      hm.dt[, Rep := tstrsplit(BioReplicate, "_",  fixed=TRUE)[3]]
      hm.dt[, xTag := do.call(paste, c(.SD, sep = " - Rep")), .SDcols=c("Condition", "Rep")]
      hm.dt[, yTag := do.call(paste, c(.SD, sep = " _ ")), .SDcols=c("Protein", "Gene")]
      
      hm.dt.w <- data.table::dcast(hm.dt, yTag ~ xTag, value.var="Abundance" )
      hm.mat <- as.matrix(hm.dt.w[,-1])
      rownames(hm.mat) <- hm.dt.w$yTag
      return(list(hm.mat,reftb))
    }
    
    
    hm.matAll <- heatmap.mat(hm.proteinAll)[[1]] # top proteins from all comparisons
    hm.matWtKo <- heatmap.mat(hm.proteinWtKo)[[1]] # top proteins all WT v. CBLB comparisons
    
    
    hm.FIX.WtKo24 <- heatmap.mat(prtFIX.WtKo24)[[1]] # fixed toplist from comparison: wt24-cblb24 (from ref. Version)
    hm.FIX.WtKo <- heatmap.mat(prtFIX.WtKo)[[1]] # fixed toplist from comparison: wt-cblb (from ref. Version)
    
    ## plot heatmaps
    hmOut.All <- hmap.plot(hm.matAll, reftb=heatmap.mat(hm.proteinAll)[[2]])
    hmOut.WtKo <- hmap.plot(hm.matWtKo, reftb=heatmap.mat(hm.proteinWtKo)[[2]])

    hmOutFIX.WtKo24 <- hmap.plot(hm.FIX.WtKo24, reftb=heatmap.mat(prtFIX.WtKo24)[[2]])
    rm <- hmap.plot(hm.FIX.WtKo, reftb=heatmap.mat(prtFIX.WtKo)[[2]])
    ################################################################# 
    
    
    
    ## save toptable and fitlists
    # if (mix3rev) {
    #   savePhrase <- paste0( "__mix3rev_",
    #                         "__FracComb_",FractComb,
    #                         "__rmNAchnnl_",allChannelsIn,"__whichPSM_",whichPSM,
    #                         "__MedPol_",medpolON,"__measurements_",meas,
    #                         "__rmMix2_",remMix2, "__onlySharedPRT_",onlySharedPRT,
    #                         "__minPptMIX_",minPptMIX, "__betweenFracNorm_", MedNorm)
    # } 

    
    ## heatmaps
    write.csv(topList, file = paste0(adir,"/loop_all/limmaTops/Top_", savePhrase, ".csv"))
    write.csv(fitList, file = paste0(adir,"/loop_all/limmaFits/Fit_", savePhrase, ".csv"))


    pdf(paste0(adir,"/loop_all/heatmap/hmapALL_",savePhrase,".pdf"), width = 16, height = 12)
    draw(hmOut.All, annotation_legend_side = "bottom")
    dev.off()

    pdf(paste0(adir,"/loop_all/heatmap/hmapWTKO_",savePhrase,".pdf"), width = 16, height = 12)
    draw(hmOut.WtKo, annotation_legend_side = "bottom")
    dev.off()
    
    
    pdf(paste0(adir,"/loop_all/heatmap/hmap_FixWTKO_",savePhrase,".pdf"), width = 16, height = 12)
    draw(hmOutFIX.WtKo, annotation_legend_side = "bottom")
    dev.off()    
    
    
    pdf(paste0(adir,"/loop_all/heatmap/hmap_FixWTKO24_",savePhrase,".pdf"), width = 16, height = 12)
    draw(hmOutFIX.WtKo24, annotation_legend_side = "bottom")
    dev.off() 
    
    
    
    ## volcano
    pdf(paste0(adir,"/loop_all/volcano/volc_all_",savePhrase,".pdf"), width = 17, height = 14)
    print(volc.all)
    dev.off()
    
    pdf(paste0(adir,"/loop_all/volcano/volc_ind_",savePhrase,".pdf"), width = 12, height = 8)
    print(volc.ind)
    dev.off()   
    
    
    ## pca
    pdf(paste0(adir,"/loop_all/pca/PCA_",savePhrase,".pdf"), width = 9, height = 7)
    print(pcaplot)
    dev.off()
    
    
    ## QC
    glist2 <- arrangeGrob(
      grobs = list(sctPlot, pDen.chnMed, pBox.chnMed, vsnTest),
      # widths = c(2,2,2,1),
      layout_matrix = rbind(c(1,1),
                            c(2,2),
                            c(3,3),
                            c(4,NA))
    )
    pdf(paste0(adir,"/loop_all/versQC/QC_",savePhrase,".pdf"), width = 7, height = 12)
    plot(glist2)
    dev.off()
    
    # glist1 <- arrangeGrob(
    #   grobs = list(MultFrac, MaxFrac, MultChnnl, perCondPlot, NAcntBar),
    #   # widths = c(2,2,2,1),
    #   layout_matrix = rbind(c(1,2),
    #                         c(3,4),
    #                         c(5,NA))
    # )
    # pdf(paste0(adir,"/loop_all/QC_plots.pdf"), width = 11, height = 14)
    # plot(glist1)
    # dev.off()
    
    
    
}}}}}}} 




