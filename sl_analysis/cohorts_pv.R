library(data.table)

setwd("/abi/data/chen")
load("cohorts.RData")
load("mut_per_cohort.RData")
load("mut1.RData")
data <- fread("/abi/data/chen/project_august/mutation.csv",data.table=F)

muts <- data.frame(pvalue)
colnames(muts) <- cohorts
muts$gene <- data$gene_name

mat = matrix(rnorm(9), 3, 3)
which(mat !=0, arr.ind = T)