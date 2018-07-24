


# rm(list=ls())

# ddir <- "/Users/farhad/_Rspace/_prj_Wolf_Isabelle/data/cblb_raw/"
# rdir <- "/Users/farhad/_Rspace/_prj_Wolf_Isabelle/results/"
# rdir <- "/Users/farhad/_Rspace/_prj_Wolf_Isabelle/scripts/"
# setwd("/Users/farhad/_Rspace/_prj_Wolf_Isabelle/")

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
library("factoextra")


# library("devtools")
# load_all(paste0("/Users/farhad/_Rspace/_loclib/MSstatsTMT_master/R"))

## ensemble mouse
ens <- useMart("ensembl", "mmusculus_gene_ensembl")




###################################
## Read in PSM data
meas <- "new"
if (meas=="new") {
  ## Reading PSM data --> new measurements. 25.06.2018 --> run 39 is missing!
  psm.raw <- as.data.table(read_excel("data/180504_P_268_IC_Pool3_7_8_SPS_180620_PSMs.xlsx"))
  ## Rename Spectrum files to match format
  psm.raw$`Spectrum File` <- gsub("_SPS","", psm.raw$`Spectrum File`)
  psm.raw$`Spectrum File` <- gsub("_tr4","", psm.raw$`Spectrum File`)
} else {
  ## Reading PSM data
  psm.raw <- as.data.table(
    read.xlsx(paste0(ddir,"171115_P_268_IC_allPools_PSMs.xlsx"),
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
psm.ano <- fread(file="data/171115_P_268_IC_Annotation.csv", sep=c(";"), header=TRUE, dec=".", key = "RAW_Filenummer")
psm.ano[, c("from", "to") := tstrsplit(RAW_Filenummer, "-",  fixed=TRUE)]
psm.ano[, c("from","to") := lapply(.SD, as.numeric), .SDcols=c("from","to")]
# psm.ano[, Run.i:= lapply(.SD, function(x) as.numeric(sub(".*\\_([^.]+)\\..*", "\\1", x))), .SDcols=Run]
psm.ano <- psm.ano[, list(Run = from:to), by = .(RAW_Filenummer, Channel, Condition, Mixture, BioReplicate)]
psm.ano[, RAW_Filenummer := NULL]
###################################



###################################
## change Charachter and Integer columns to Factor
## PSM table:
psm.dt <- char2fact(psm.dt)
## Anootation table:
psm.ano <- char2fact(psm.ano)
psm.ano <- int2fact(psm.ano)
###################################


