### R code from vignette source 'RTopper.Rnw'

###################################################
### code chunk number 1: start
###################################################
options(width=50)
rm(list=ls())


###################################################
### code chunk number 2: loadData
###################################################
library(RTopper)
data(exampleData)
ls()
class(dat)
names(dat)
sapply(dat,class)
sapply(dat,dim)
dim(pheno)
str(pheno)


###################################################
### code chunk number 3: strData
###################################################
###data structure
lapply(dat,function(x) head(x)[,1:3])
sum(rownames(dat[[1]])%in%rownames(dat[[2]]))
sum(rownames(dat[[2]])%in%rownames(dat[[3]]))


###################################################
### code chunk number 4: metaData
###################################################
library(org.Hs.eg.db)
org.Hs.eg()


###################################################
### code chunk number 5: listFGS
###################################################
kegg <- as.list(org.Hs.egPATH2EG)
go <- as.list(org.Hs.egGO2ALLEGS)
length(kegg)
length(go)
str(kegg[1:5])
names(kegg)[1:5]
str(go[1:5])
names(go)[1:5]


###################################################
### code chunk number 6: convertIDs
###################################################
kegg <- lapply(kegg[sample(1:length(kegg),5)],function(x) unique(unlist(mget(x,org.Hs.egSYMBOL))))
go <- lapply(go[sample(1:length(go),5)],function(x) unique(unlist(mget(x,org.Hs.egSYMBOL))))
str(kegg)
str(go)


###################################################
### code chunk number 7: annotateFGS
###################################################
library(KEGG.db)
KEGG()
names(kegg) <- paste(names(kegg),unlist(mget(names(kegg),KEGGPATHID2NAME)),sep=".")
names(kegg)


###################################################
### code chunk number 8: listFGS
###################################################
library(GO.db)
GO()
names(go) <- paste(names(go),Term(names(go)),sep=".")
names(go)


###################################################
### code chunk number 9: listFGS
###################################################
fgsList <- list(go=go,kegg=kegg)


###################################################
### code chunk number 10: convertToDr
###################################################
dataDr <- convertToDr(dat, pheno, 4)
class(dataDr)
length(dataDr)
names(dataDr)[1:5]
str(dataDr[1:2])


###################################################
### code chunk number 11: integratedScore
###################################################
bicStatInt <- computeDrStat(dataDr, columns = c(1:4), method="bic", integrate = TRUE)
names(bicStatInt)
str(bicStatInt)


###################################################
### code chunk number 12: separateScore
###################################################
bicStatSep <- computeDrStat(dataDr, columns = c(1:4), method="bic", integrate = FALSE)
names(bicStatSep)
str(bicStatSep)


###################################################
### code chunk number 13: runGSEbatchArgs
###################################################
args(runBatchGSE)


###################################################
### code chunk number 14: runBatchGSE.int1
###################################################
gseABS.int <- runBatchGSE(dataList=bicStatInt, fgsList=fgsList)
gseABS.int <- runBatchGSE(dataList=bicStatInt, fgsList=fgsList,
                          absolute=TRUE, type="f", alternative="mixed")


###################################################
### code chunk number 15: runBatchGSE.int2
###################################################
gseUP.int <- runBatchGSE(dataList=bicStatInt, fgsList=fgsList,
                         absolute=FALSE, type="t", alternative="up")
gseDW.int <- runBatchGSE(dataList=bicStatInt, fgsList=fgsList,
                         absolute=FALSE, type="t", alternative="down")
gseBOTH.int <- runBatchGSE(dataList=bicStatInt, fgsList=fgsList,
                           absolute=FALSE, type="t", alternative="either")


###################################################
### code chunk number 16: runBatchGSE.int3
###################################################
gseABSsim.int <- runBatchGSE(dataList=bicStatInt, fgsList=fgsList,
                             absolute=TRUE, type="f", alternative="mixed",
                             ranks.only=FALSE, nsim=1000)
gseUPsim.int <- runBatchGSE(dataList=bicStatInt, fgsList=fgsList,
                            absolute=FALSE, type="t", alternative="up",
                            ranks.only=FALSE, nsim=1000)


###################################################
### code chunk number 17: runBatchGSE.format1
###################################################
str(gseUP.int)
gseABSsim.int


###################################################
### code chunk number 18: runBatchGSE.altFunc
###################################################
library(limma)
gseUP.int.2 <- runBatchGSE(dataList=bicStatInt, fgsList=fgsList,
                           absolute=FALSE, gseFunc=wilcoxGST, alternative="up")


###################################################
### code chunk number 19: runBatchGSE.format2
###################################################
str(gseUP.int.2)
all(gseUP.int.2$go==gseUP.int$go)


###################################################
### code chunk number 20: runBatchGSE.altFunc2
###################################################
gseFunc <- function (selected, statistics, threshold) {
  diffExpGenes <- statistics > threshold
  tab <- table(diffExpGenes, selected)
  pVal <- fisher.test(tab)[["p.value"]]
}
gseUP.int.3 <- runBatchGSE(dataList=bicStatInt, fgsList=fgsList,
                           absolute=FALSE, gseFunc=gseFunc, threshold=7.5)


###################################################
### code chunk number 21: runBatchGSE.format3
###################################################
str(gseUP.int.3)
data.frame(fisher=gseUP.int.3$integrated$kegg,wilcoxon=gseUP.int$integrated$kegg)


###################################################
### code chunk number 22: runBatchGSE.separate
###################################################
gseABS.sep <- runBatchGSE(dataList=bicStatSep, fgsList=fgsList)


###################################################
### code chunk number 23: combineGSE
###################################################
gseABS.geoMean.sep <- combineGSE(gseABS.sep, method="geometricMean")
gseABS.max.sep <- combineGSE(gseABS.sep, method="max")


###################################################
### code chunk number 24: combineGSE.format
###################################################
names(gseABS.sep)
str(gseABS.sep)
str(gseABS.geoMean.sep)
gseABS.geoMean.sep


###################################################
### code chunk number 25: adjustP
###################################################
gseABS.int.BH <- adjustPvalGSE(gseABS.int)
gseABS.int.holm <- adjustPvalGSE(gseABS.int, proc = "Holm")


###################################################
### code chunk number 26: adjusted.format
###################################################
names(gseABS.int.BH)
names(gseABS.int.holm)
str(gseABS.int.BH)
str(gseABS.int.holm)


###################################################
### code chunk number 27: sessioInfo
###################################################
sessionInfo()

