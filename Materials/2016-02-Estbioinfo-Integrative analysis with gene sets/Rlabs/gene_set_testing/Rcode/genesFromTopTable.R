#' Extract a list (a character vector indeed) of genes from a topTable outputted by the limma software
#' @param aTopTable
#' @param filename
#' @param entrezOnly
#' @param uniqueIds
#' @param adjOrrawP
#' @param Pcutoff
#' @param FCcutoff
#' @param id2select
#' @param cols2select
#' @keywords
#' @seealso
#' @export
#' @examples
#'

genesFromTopTable <- function (aTopTable,
                               filename=NULL, 
                               entrezOnly = TRUE, 
                               uniqueIds=TRUE, 
                               adjOrrawP = "adj",
                               Pcutoff = 0.05,
                               FCcutoff = 1,
                               id2Select = "ENTREZ",  # c("ENTREZ","SYMBOL", "other") 
                                                      # Alias ("Entrez", "EntrezsA", "Symbols", "SymbolsA") are accepted
                               cols2Select = 2){
  
  if (! is.null(filename)){
    topTab <- read.csv(file.path(listsDir, filename), head=TRUE, sep=";", dec=",", row.names=1)
  }else{
    topTab=aTopTable
  }
  dim(topTab)
  colnames(topTab)
  if (entrezOnly) {
    selectedEntrez <- !is.na(topTab[,"EntrezsA"])
    topTab <- topTab[selectedEntrez,]
  }
  dim(topTab)
  if (Pcutoff < 1){
    if (adjOrrawP=="adj"){
      selectedP <- topTab[,"adj.P.Val"] < Pcutoff
    }else{
      selectedP <- topTab[,"P.Value"] < Pcutoff
    }
    topTab<- topTab[selectedP, ]
  }
  dim(topTab)
  if (FCcutoff > 0){
    selectedFC <-(abs(topTab[,"logFC"]) > FCcutoff)
    topTab<- topTab[selectedFC, ]
    dim(topTab)
  }
  
  if (toupper(substr(id2Select,1,6))=="SYMBOL"){
    geneList <- topTab[,"SymbolsA"]
  }else{
    if (toupper(substr(id2Select,1,6))=="ENTREZ"){
      geneList <- topTab[,"EntrezsA"]
    }else{
      if(cols2Select> 0){
        geneList <- topTab[,cols2Select]
      }else{
        geneList <-rownames(topTab)
      }}}
  # length(geneList)
  if(uniqueIds) geneList <- unique(geneList)
  return(geneList)
}


