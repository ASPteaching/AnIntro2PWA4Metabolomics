#Will use RTopper and limma packages to perform single platform and integrated gene sets analysis. 
if (!(require(RTopper))){
  source("http://bioconductor.org/biocLite.R")
  biocLite("RTopper")
}
require(RTopper)

##Load sample TCGA data that comes with the package
data(exampleData)

## view loaded objects
objects()

## or
ls()

## you should find objects dat and pheno
## Veiw what is inside the objects:
names(dat)
str(dat)
str(pheno)

sapply(dat,dim)

##access Agilent-based gene expression data and print several rows and columns
dat$dat.affy[1:5, 1:5]
##access affymetrix-based gene expression data and print several rows and columns
dat$dat.agilent[1:5, 1:5]

### WE SEE: SAME COLUMN NAMES, SAME ROW NAMES, but DIFFERENT EXPRESSION VALUES

##acces first several rows of the phenotype data

pheno[1:10,]

## Perform GSA for gene expression.

#first let's compute gene to phenotype association scores -- we will use an absolute value of the moderated t-statistic prodcued by package limma
library(limma)
#create a design matrix: 
design.matrix = model.matrix(~pheno$Class)

#fit linear model and extract moderated t
x = lmFit(dat$dat.affy, design.matrix)
fit = eBayes(x)
affy.de = topTable(fit, coef=2, n=nrow(fit$t), adjust.method="BH") 
affy.de[1:10,]

#load R-object with KEGG gene sets
load("Msig_c2_cp_kegg_v3.rda")
# stud loaded object and view first 5 sets from a loaded object
str(Msig_c2_cp_kegg_v3)
Msig_c2_cp_kegg_v3$data[1:5]

## perform gene sets testing using limma's geneSetTest function 
## view help page for the function 
?geneSetTest

## we will run test for KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION pathway
## compare
Msig_c2_cp_kegg_v3$data[89]
#and
Msig_c2_cp_kegg_v3$data[[89]]

set = Msig_c2_cp_kegg_v3$data[[89]]

## create indicator for genes in the set: 

index = rownames(affy.de)%in%set

geneSetTest(index, abs(affy.de$t))

#Is our result surprising?
plot(density(abs(affy.de$t)[!index]), col="black", lwd=2, xlim=c(-0.5,5))
lines(density(abs(affy.de$t)[index]), col="red", lwd=2)

########### Example of integrative GSA

#Convert data to a format required to compute integrated association scores
dataDr <- convertToDr(dat, pheno, 4)

#Compute individual scores using logistic regression (for a single platform it should give similar rankings, as obtained by using moderated t from limma)

devStatSep <- computeDrStat(dataDr, columns = c(1:4), method="dev", integrate = FALSE)
names(devStatSep)

devStatInt <- computeDrStat(dataDr, columns = c(1:4), method="dev", integrate = TRUE)
names(devStatInt)

#Run Gene sets analysis for each individual data type
gseABS.Sep<- runBatchGSE(dataList=devStatSep, fgsList=list(KEGG_pathways=Msig_c2_cp_kegg_v3$data),
                         absolute=TRUE, type="f", alternative="mixed")
str(gseABS.Sep)

#Run Gene sets analysis for integrated score
gseABS.int <- runBatchGSE(dataList=devStatInt, fgsList=list(KEGG_pathways=Msig_c2_cp_kegg_v3$data),
                          absolute=TRUE, type="f", alternative="mixed")
#Multiple testing correction
gseABS.sep.BH <- adjustPvalGSE(gseABS.Sep)
gseABS.int.BH <- adjustPvalGSE(gseABS.int)

