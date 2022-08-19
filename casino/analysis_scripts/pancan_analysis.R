library("ggplot2")
library("reshape2")
library("plyr")
library("scales")

options(echo=TRUE)
args<-commandArgs(trailingOnly = TRUE)

n = 30

#cohs <- c("PBCA-DE","PACA-AU","PACA-CA","PAEN-AU","SARC-US","OV-AU","OV-US","BRCA-EU",
#          "BRCA-UK","PRAD-UK","PRAD-US","EOPC-DE")

cohort <- args[1]
#cohort <- ("PRAD-UK,PRAD-US,EOPC-DE")
cohs <- strsplit(cohort, split=",")
cohs <- unlist(cohs)

elements_list = list()

for (i in 1:length(cohs)){
  
  #print(cohs[i])
  score_file <- read.table(paste("/ibios/co02/chen/new_casino/reports/enhancer_",cohs[i],".txt",sep=""),header=T)
  top_n <- head(score_file[order(-score_file$snv_rs),],n)
  #gene <- list(as.character(top_n$gene))
  #gene_list <- c(gene_list,gene)
  rows <- row.names(top_n)
  elements_list <- c(elements_list,rows)
}


elements <- unique(unlist(elements_list))

#for (i in 1:length(head(genes))){
#  if (head(genes)[i] %in% genes){  
#    print(i)
#    print(head(genes)[i])
    
#  }
  
#}

score_file <- read.table(paste("/ibios/co02/chen/new_casino/reports/enhancer_",cohs[i],".txt",sep=""),header=T)
selected <- score_file[which(row.names(score_file) %in% elements),]
target <- data.frame(selected[,1:4])

for (i in 1:length(cohs)){
  
  score_file <- read.table(paste("/ibios/co02/chen/new_casino/reports/enhancer_",cohs[i],".txt",sep=""),header=T)
  selected <- score_file[which(row.names(score_file) %in% elements),]
  target[,eval(cohs[i])]<-selected$snv_rs
}

target.m <- target[,-c(1,2,3)]
target.s <- melt(target.m)
target.s <- ddply(target.s, .(variable), transform, rescale = rescale(value))
(p <- ggplot(target.s, aes(variable, gene)) + geom_tile(aes(fill = rescale),colour = "white") 
 + scale_fill_gradient(low = "white", high = "steelblue"))
#p + theme(axis.text.y = element_blank())

jpeg(filename=paste("/home/hongc/pancan/snvs/enhancer2_",cohort,".jpeg",sep=""))
p + theme(axis.text.y = element_blank())
dev.off()
