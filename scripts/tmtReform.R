



###################################
## Read in PSM data
meas <- "new"
if (meas=="new") {
  ## Reading PSM data --> new measurements. 25.06.2018 --> Fraction 39 is missing!
  psm.raw <- as.data.table(read_excel(paste0(ddir,"180504_P_268_IC_Pool3_7_8_SPS_180620_PSMs.xlsx")))
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
psm.tmp[, Fraction:= lapply(.SD, function(x) as.integer(sub(".*\\_([^.]+)\\..*", "\\1", x))), .SDcols="Spectrum.File"]
i.x <- which(names(psm.tmp) %like% "Abundance")
setnames(psm.tmp, old = names(psm.tmp)[i.x], new = gsub("Abundance:.", "", names(psm.tmp)[i.x]))
psm.dt <- melt.data.table(psm.tmp, id.vars = names(psm.tmp)[!names(psm.tmp) %in% c("126", "127", "128", "129", "130", "131")],
                          variable.name="Channel", 
                          value.name="Intensity")
###################################



###################################
## Reading GENE IDS 
dt.ID <- as.data.table(read_excel(
  paste0(ddir,"171115_P_268_IC_replicates3_7_8_171129_IDs.xlsx"), 
  sheet = 2), key="Accession")
dt.ID <- dt.ID[,list(Accession, `Gene Symbol`)]
setnames(dt.ID, old = "Gene Symbol", new = "Gene")




###################################



###################################
## Reading PSM Annotation table and reformatting 
psm.ano <- fread(file=paste0(ddir,"171115_P_268_IC_Annotation.csv"), sep=c(";"), 
                 header=TRUE, dec=".", key = "RAW_Filenummer")
psm.ano[, c("from", "to") := tstrsplit(RAW_Filenummer, "-",  fixed=TRUE)]
psm.ano[, c("from","to") := lapply(.SD, as.numeric), .SDcols=c("from","to")]
# psm.ano[, Fraction.i:= lapply(.SD, function(x) as.numeric(sub(".*\\_([^.]+)\\..*", "\\1", x))), .SDcols=Fraction]
psm.ano <- psm.ano[, list(Fraction = from:to), by = .(RAW_Filenummer, Channel, Condition, Mixture, BioReplicate)]
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



########################################################
################################### Filterring procedure
## Merge with annotation
psm.ready <- merge(psm.dt, psm.ano, by = c("Fraction","Channel"), allow.cartesian=TRUE)
psm.ready <-  psm.ready[, c("Fraction", "Channel", "Master.Protein.Accessions", "#.Protein.Groups", "Annotated.Sequence", "Charge", 
                    "Ions.Score", "Rank", "MS.Order", "File.ID", "RT.[min]", "Condition", "Quan.Info",
                    "Mixture", "BioReplicate", "Intensity")]
setnames(psm.ready, old = c("Master.Protein.Accessions", "Annotated.Sequence", "#.Protein.Groups"), 
         new = c("Protein", "Peptide", "numProtein") )
psm.ready[Intensity < 1, Intensity := 1]

####### STEP 1: Removing shared proteins
## Number of protein accessions in Protein --> this should be equal to 1.
psm.ready[, cntProt := lengths(strsplit(as.character(Protein), ";"))]
# Remove rows with Protein == "sp"
psm.ready <- psm.ready[!Protein=="sp", ] 
# Remove rows with more than 1 protein accession and no "sp" (shared protein)
psm.ready <- psm.ready[cntProt>1 & Protein %like% "sp" | cntProt==1, ]
## For Protein containing "sp" and one protein accession, remove "sp" and stay with protein accession.
## If the number of Protein accessions after removing "sp" is more than 1 --> shared protein --> remove row. 
psm.ready[,Protein := gsub("sp;","", Protein, fixed = TRUE)]
psm.ready[,Protein := gsub("; sp","", Protein, fixed = TRUE)]
psm.ready[,Protein := gsub(" ","", Protein, fixed = TRUE)]
psm.ready <- psm.ready[!Protein=="sp", ] 

psm.ready[, cntProt := lengths(strsplit(as.character(Protein), ";"))]
psm.ready <- psm.ready[cntProt==1]
## The process remove 150462 rows --> 11% of the data
# psm.ready[,Quan.Info:=factor("Unique")] #Now Quant.info can be marked as unique. ## this id needed for MSstatsTMT

psm.ready[,cntProt:=NULL]
## define Feature and Proptide
psm.ready[, Feature:= do.call(paste, c(.SD, sep = "_")), .SDcols=c("Protein", "Peptide", "Charge")]
psm.ready[, Proptide:= do.call(paste, c(.SD, sep = "_")), .SDcols=c("Protein", "Peptide")]
## Define Abundance as log2(Intensity)
psm.ready[, Abundance0 := log2(Intensity)]
# Continue with Abundance: Andreas 6.6.2018






################################################################# 
## extract bioloical parameters  ..> genes
## ensemble mouse
ens <- useMart("ensembl", "mmusculus_gene_ensembl", host="www.ensembl.org" )
# mmusculus <- useMart(host="useast.ensembl.org", 
#                      biomart="ENSEMBL_MART_ENSEMBL", 
#                      dataset="mmusculus_gene_ensembl")

# uniProt <- useMart("unimart", dataset="uniprot")


trans = getBM(attributes = c("uniprotswissprot", "mgi_symbol", "ensembl_gene_id", "mgi_symbol",
                             "external_gene_name", "description"), 
              # filters = "uniprotswissprot", values = "Q8VCN9",
              filters = "uniprotswissprot", values = unique(psm.ready$Protein),
              uniqueRows=FALSE, mart = ens, verbose = FALSE)


trans <- as.data.table(trans)

setkey(psm.ready, "Protein")
setkey(trans, "uniprotswissprot")
psm.ready[trans, Gene.ens := external_gene_name]
psm.ready[trans, Gene.mgi := mgi_symbol]

setkey(dt.ID, "Accession")
psm.ready[dt.ID, Gene.id := Gene]

psm.ready$Genetmp <- ifelse(!is.na(psm.ready$Gene.id), psm.ready$Gene.id, psm.ready$Gene.ens)
psm.ready$Gene <- ifelse(!is.na(psm.ready$Genetmp), psm.ready$Genetmp, psm.ready$Gene.mgi)


################################################################# 
## add run
  tmp <- unique(psm.ready[,list(Fraction, Mixture)])
  tmp <- tmp[order(Mixture, Fraction)]
  tmp[ , Frac.i := seq_along(Fraction), by=Mixture]
  tmp[ , Run := paste0(seq_along(Fraction), ".", Mixture), by=Mixture]
  
  psm.ready <- merge(psm.ready, tmp, by=c("Fraction","Mixture"), all.x = TRUE)
  

 # RawBox <- 
 #  ggplot(psm.ready, aes(x=Fraction, y=Abundance0, col=Condition)) + 
 #    geom_boxplot() +
 #    scale_color_brewer(palette = "Dark2") + 
 #    theme_bw() + 
 #    # geom_bar(aes(y = (..count..)/sum(..count..)))
 #    # geom_bar(stat = "count") +
 #    facet_wrap(~Mixture, scales = "free_x") 
 #  
 #  
 # RawBar <-
 #   ggplot(psm.ready, aes(x=Fraction, col=Condition)) +
 #   scale_color_brewer(palette = "Dark2") +
 #   theme_bw() +
 #   # geom_bar(aes(y = (..count..)/sum(..count..)))
 #   geom_bar(stat = "count", position = "dodge2") +
 #   facet_wrap(~Mixture, scales = "free_x")
 


  


