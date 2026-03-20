library(MetaboAnalystR)
mSet <- InitDataObjects("list", "pathora", FALSE, default.dpi=300)
cmpd.vec<-c("Acetoacetic acid","Beta-Alanine","Creatine","Dimethylglycine","Fumaric acid","Glycine","Homocysteine","L-Cysteine","L-Isolucine","L-Phenylalanine","L-Serine","L-Threonine","L-Tyrosine","L-Valine","Phenylpyruvic acid","Propionic acid","Pyruvic acid","Sarcosine","Arsenic","Benzene","Caffeic acid","Cotinine","Cadmium","Lead","Thiocyanate")
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "name");
mSet<-CreateMappingResultTable(mSet)
mSet<-CreateMappingResultTable(mSet)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 2);
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_0_", "net", "png", 150, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_0_", "png", 150, width=NA)
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_1_", "net", "png", 150, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_1_", "png", 150, width=NA)
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_2_", "net", "png", 150, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_2_", "png", 150, width=NA)
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_3_", "net", "png", 150, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_3_", "png", 150, width=NA)
