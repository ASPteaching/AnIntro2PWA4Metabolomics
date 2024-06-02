# Cargamos los datos preprocesados para el GSA

load(file="data4GSEA.Rda")

# Preparamos la matriz binaria para los GeneSets

GSMatrix <- matrix(NA, nrow=length(entrezs), ncol=length(gsEntrez))
for (i in 1:length(gsEntrez)){
    GSMatrix[,i] <- entrezs %in% gsEntrez[[i]] # binaryGenSets[[i]]
}
GSMatrix <- GSMatrix*1
colnames(GSMatrix)<- names(gsEntrez)
rownames(GSMatrix)<- entrezs
head(GSMatrix)

# Cada columna de esta matriz tiene unos en las filas correspondientes a sus genes
# Comprobamos que da lo mismo que en los Genesets individuales

require(naturalsort)
for (i in 1:ncol(GSMatrix)){
    GSm<- naturalsort(rownames(GSMatrix[GSMatrix[,i]!=0,]))
    gsE<-naturalsort(gsEntrez[[i]])
    sumEq <- sum(GSm!=gsE)
    cat(sumEq, " differences found in Gene Set ", i, " ", colnames(GSMatrix)[i], "\n")
}

    
# Para empezar calculamos un test t para cada gen
  
require(genefilter)
rownames(expresHvsL)<-entrezs
tScores <-rowttests(expresHvsL, groupsHvsL, tstatOnly=TRUE)
head(tScores)
    
# Para decidir si un GeneSet esta sobre(sub) representado haremos un test de Wilcoxon sobre los resultados de los ttest
  
for (i in (1:ncol(GSMatrix))){
    text2show <- paste("Gene Set ", i, ": ", colnames(GSMatrix)[i])
    cat(text2show,"\n")
    cat(paste(rep("=", rep(nchar(text2show))), collapse=""),"\n")
    pOver <-wilcox.test(tScores$statistic~GSMatrix[,i], alternative="greater")$p.value
    pUnder <-wilcox.test(tScores$statistic~GSMatrix[,i], alternative="less")$p.value
    cat("\tTest for OverExpression :", pOver,"\n")
    cat("\tTest for UnderExpression :", pUnder,"\n")
}

# Los resultados no difieren mucho de los que se obtiene con un test "estandar" de GSEA
con    <- textConnection('stdout', 'wr', local = TRUE)
sink(con)
GSA.obj<-GSA(expresHvsL, groupsHvsL, genenames=entrezs, genesets=gsEntrez,  
             resp.type="Two class unpaired", nperms=1000)  # Provar amb 1000?
sink()
text2show <- paste("Conjuntos de genes que aparecen en la comparacion",                      "H vs L")
cat("\n",text2show,"\n")
cat(paste(rep("=", rep(nchar(text2show))), collapse=""),"\n")
genesetsNames <- names(gsEntrez)
show(GSA.listsets(GSA.obj, geneset.names=genesetsNames, FDRcut=.5))
