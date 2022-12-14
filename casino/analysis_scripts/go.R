source("http://www.bioconductor.org/biocLite.R")
biocLite(c("biomaRt", "topGO", "org.Mm.eg.db"))
library(biomaRt)
library(org.Mm.eg.db)
library(topGO) 

all_genes <- rep(0, each=145)

names(all_genes) <- c("ENSMUSG00000031167", "ENSMUSG00000022443", "ENSMUSG00000055762", "ENSMUSG00000031328", 
                      "ENSMUSG00000026728", "ENSMUSG00000029722", "ENSMUSG00000059208", "ENSMUSG00000017404", 
                      "ENSMUSG00000003970", "ENSMUSG00000019210", "ENSMUSG00000024966", "ENSMUSG00000024197", 
                      "ENSMUSG00000030847", "ENSMUSG00000004264", "ENSMUSG00000041959", "ENSMUSG00000032643", 
                      "ENSMUSG00000042520", "ENSMUSG00000028484", "ENSMUSG00000022186", "ENSMUSG00000027593", 
                      "ENSMUSG00000019795", "ENSMUSG00000021134", "ENSMUSG00000026956", "ENSMUSG00000022982", 
                      "ENSMUSG00000027162", "ENSMUSG00000024668", "ENSMUSG00000000168", "ENSMUSG00000022234", 
                      "ENSMUSG00000024991", "ENSMUSG00000018697")


all_genes <- c(1,1)
names(all_genes) <- c("TP53","ERG")
  
GOdata <- new("topGOdata", ontology = "BP", allGenes = all_genes, geneSel = function(p) p < 1e-2, 
              description = "Test", annot = annFUN.org, mapping="org.Mm.eg.db", ID="Ensembl")

#GOdata <- new("topGOdata", ontology = "BP", allGenes = all_genes, geneSel = function(p) p < 1e-2, 
              description = "Test", annot = annFUN.org, mapping="org.Mm.eg.db", ID="GeneName")

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
GenTable(GOdata, classicFisher = resultFisher, topNodes = 10)
