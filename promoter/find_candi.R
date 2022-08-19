library(data.table)
library(ggplot2)


data <- fread(file="/ibios/co02/chen/pcawg_2018/result/matrix/promoter/select/pm_snv.csv",data.table=F)

cood <- fread(file="/ibios/co02/chen/pcawg_2018/result/matrix/promoter/BLCA-US2.bed",data.table=F)
cood <- cood[,c(1:4)]

cood_ensg <- fread(file="/home/hongc/PhD/promoters/code/promoter_match.txt",data.table=F,header = F)

datasum <- rowSums(data)
cood$datasum <- datasum
colnames(cood_ensg) <- c("ensg","hg")
cood$ensg <- cood_ensg$ensg
mut <- cood[order(-cood$datasum),]
write.table(mut,file="/ibios/co02/chen/pcawg_2018/result/matrix/promoter/report/stats.txt",quote=F,row.names=F,col.names=T,sep="\t")

full.table <- cbind(cood,data)
write.table(full.table,file="/ibios/co02/chen/pcawg_2018/result/matrix/promoter/report/full_prm_snv.txt",quote=F,row.names=F,col.names=T,sep="\t")


## test ##
ds <- cbind(cood,data)
df <- ds[(ds$gene=="AQP1A"|ds$gene=="AQP1B"),]
df <- df[,5:ncol(df)]
df[3,] <- colSums(df)
df2 <- data.frame(t(df))
colnames(df2) <- c("a1","a2","sum")
df2 <- df2[df2$sum==1,]
