metanr_packages <- function(){
  
  metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest","RBGL","edgeR","fgsea","devtools","crmn","httr","qs")
  
  list_installed <- installed.packages()
  
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  
  if(length(new_pkgs)!=0){
    
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(new_pkgs)
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}

metanr_packages()

library(MetaboAnalystR)
mSet <- InitDataObjects("list", "pathora", FALSE) default.dpi=300)
cmpd.vec <- c("Acetoacetic acid","Beta-Alanine","Creatine","Dimethylglycine","Fumaric acid",
              "Glycine","Homocysteine","L-Cysteine","L-Isolucine","L-Phenylalanine","L-Serine","L-Threonine","L-Tyrosine","L-Valine","Phenylpyruvic acid","Propionic acid", "Pyruvic acid","Sarcosine","Arsenic","Benzene","Caffeic acid","Cotinine", "Cadmium","Lead","Thiocyanate")

mSet <- Setup.MapData(mSet, cmpd.vec)
mSet <- CrossReferencing(mSet, "name")
sessionInfo()