

ddir <- "/Users/farhad/_Rspace/cblb_wtVkn/data/"
rdir <- "/Users/farhad/_Rspace/cblb_wtVkn/results/"

library("data.table")
library("dplyr")
library("ggplot2")
library("corrplot")
library(RColorBrewer)

# library("devtools")
# load_all(paste0("/Users/farhad/_Rspace/_loclib/MSstatsTMT_master/R"))


inp <- readRDS( file = "required_input.rds")
qnt <- readRDS( file = "quant_byprotein_NORM.rds")

qnt.w <- dcast(qnt, Protein ~ Condition + BioReplicate, value.var="Abundance.Norm")
qnt.w[, nNA := Reduce(`+`, lapply(.SD,function(x) is.na(x))), .SDcols=names(qnt.w)[-c(1,2)]]


qnt.m <- as.matrix(qnt.w[,-c("Protein","nNA")])

M <- cor(qnt.m, method = "pearson", use = "pairwise.complete.obs")
res <- cor.mtest(M, conf.level = .95)

corrplot(M, p.mat = res$p, method = "color", type = "lower",
         sig.level = c(.001, .01, .05), pch.cex = 1, pch.col = "red",
         tl.col="black", tl.cex = 0.9, tl.srt = 45,
         order = "hclust", cl.lim = c(0,1),
         number.cex = .8, addCoef.col = "white", diag = FALSE)

