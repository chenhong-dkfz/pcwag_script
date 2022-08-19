library(cometExactTest)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
gene.name <- args[1]
muts <- args[2]

setwd("/abi/data")
data <- fread(paste0("/abi/data/chen/project_august/mutation25/",gene.name,"/",muts,".txt"),data.table=F)

# file input : n_total, n_mut1, n_mut2, n_intersect
# comet input, c(n_total,n_mut1, n_mut2, n_intersect)
# phyper input, n_intersect, n_mut1, n_total-n_mut1, n_mut2


hyperdis <- function(a){
  total <- a[[1]]
  mut1 <- a[[2]]
  mut2 <- a[[3]]
  inter <- a[[4]]
  ndel <- total - mut1
  p.value <- phyper(inter,mut1,ndel,mut2,lower.tail = T)
}

# 
# comdis <- function(a){
#   total <- a[[1]]
#   mut1 <- a[[2]]
#   mut2 <- a[[3]]
#   inter <- a[[4]]
#   ndel <- total - mut1
#   p.value <- comet_exact_test(c(total,mut1,mut2,inter))
# }
# 
# 


mut0 <- apply(data, 1, hyperdis)
#mut1 <- apply(data2, 1, comdis)
save(mut0,file=paste0("/abi/data/chen/project_august/mutation25/result/",gene.name,"/",muts,".Rdata"))