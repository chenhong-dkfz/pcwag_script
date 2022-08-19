library(data.table)
library(stringr)


setwd("/abi/data/chen/promoter_project/")

pm <- fread("fout.txt",sep=",",data.table = F)
coding <- fread("fout2.txt",data.table = F,header = F)

pm$V3 <- 0
pm[pm$V1 %in% coding$V1,3] <- 1

ts <- str_split_fixed(pm$V1, "_", 2)
pm <- cbind(pm,ts)
names(pm) <- c("combi","zscore","amplification","pid","gene")

cgc.table <- fread("/home/hongc/PhD/basic/cgc.txt",data.table = F)
cohort.data <- fread("/abi/data/chen/project_august/pcawg_ta_project.txt",data.table = F)
names(cgc.table)[1] <- "gene"
names(cohort.data)[1] <- "pid"

pm$cgc <- 0
pm[pm$gene %in% cgc.table$gene,6] <- 1
pm2 <- merge(pm,cohort.data,by="pid")
pm <- pm2
pm <- pm[,-2]
pm.cgc <- pm[pm$cgc==1,]
pm2 <- merge(pm.cgc,cgc.table,by="gene")
pm.cgc <- pm2
write.table(pm,file="promoter_overexpression.txt",sep="\t",quote = F,row.names = F,col.names = T)
write.table(pm.cgc,file="cgc_promoter_overexpression.txt",sep="\t",quote = F,row.names = F,col.names = T)




