maxVarByID <- function (x, genes, numericCols=NULL){
  lincs<-cbind(P=rownames(x), G=genes)  
  geneNames <- unique(genes)
  numGenes<-length(geneNames)
  maxVars<-matrix(0, nrow=numGenes, ncol=ncol(x))
  rownames(maxVars) <- geneNames
  colnames(maxVars) <- colnames(x)
  i<-0
  for(g in geneNames){
    i<-i+1
    subsX <- as.matrix(x[lincs[,2]==g,])
    if (nrow(subsX)==1){
      maxVars[i,]<- subsX
    }else{
      sds <- apply(subsX, 1, sd)
      maxIdx <- which(sds==max(sds))
      maxVars[i,] <-subsX[maxIdx,]
    }
  }
  return(maxVars)
}

# test
# x<- matrix(c(1,2,2,2,1, 1,3,5,3,1, 1,1,2,1,1, 1,5,8,5,1, 1,1,2,2,1), nrow=5, byrow=TRUE)
# rownames(x)<-paste("p",1:5, sep="");x
# genes <- paste("g", c(1,2,2,3,3), sep="")
# aggregate(x, by=list(genes), FUN=mean)
# maxVarByID (x, genes)
