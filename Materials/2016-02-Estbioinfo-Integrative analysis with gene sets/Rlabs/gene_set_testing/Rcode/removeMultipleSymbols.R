# Gene Symbols are often provided with all their synonims, separated by a double bacslash '//'
# Informative as this can be this makes lists unmanageable for setting data rownames, printing or for input to other analyses
# An alternative is to keep only the first symbol with the hope that this is the most used.

removeMultipleIDs <- function (multipleIDsSymbol="//", myAnotTable){
  if (length(grep(multipleIDsSymbol, myAnotTable))>0) {
    # if (runMulticore ==1 || runMulticore ==3) { 
    # myAnotTable <- unlist(mclapply(strsplit(myAnotTable, multipleIDsSymbol),
    #                                 my.fun <- function(x){return(x[1])}))
    #} else {
      myAnotTable <- unlist(lapply(strsplit(myAnotTable, multipleIDsSymbol),
                                   my.fun <- function(x){return(x[1])}))         
    #}
  }
}

# Example
# symbs <-as.character(read.table(file="symbolsList.txt", head=FALSE)[,1])
# cleanSymbs <- removeMultipleIDs (myAnotTable=symbs)
# data.frame(symbs, cleanSymbs)
