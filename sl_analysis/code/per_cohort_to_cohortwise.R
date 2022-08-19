# this script convert pasted 4 values to one

library(data.table)
library(ComplexHeatmap)
library(dplyr)
library(easyGgplot2)
library(grid)
library(gridExtra)
library(cometExactTest)
library(data.table)


hyperdis <- function(a){
  total <- a[[1]]
  mut1 <- a[[2]]
  mut2 <- a[[3]]
  inter <- a[[4]]
  ndel <- total - mut1
  p.value <- phyper(inter,mut1,ndel,mut2,lower.tail = T)
  return(p.value)
}

fishertest <- function(a){
  inter <- a[[4]]
  mut1s <- a[[2]]-a[[4]]
  mut2s <- a[[3]]-a[[4]]
  nonmut <- a[[1]]-a[[2]]-a[[3]]+a[[4]]
  mtx <- matrix(c(inter,mut1s,mut2s,nonmut),nrow=2)
  res <- fisher.test(mtx)
  p.value <- res$p.value
  odds <- res$estimate[[1]]
  return(c(p.value,odds))
}


df.mut <- fread("/abi/data/chen/project_august/mutation25_tial1/mut_mut.txt",data.table = F)

df.1 <- df.mut[,seq(1, ncol(df.mut), 4) ]
df.2 <- df.mut[,seq(2, ncol(df.mut), 4) ]
df.3 <- df.mut[,seq(3, ncol(df.mut), 4) ]
df.4 <- df.mut[,seq(4, ncol(df.mut), 4) ]

vt1 <- rowSums(df.1)
vt2 <- rowSums(df.2)
vt3 <- rowSums(df.3)
vt4 <- rowSums(df.4)

df.merge <- data.frame(cbind(vt1,vt2,vt3,vt4))
data <- df.merge
pvalue <- apply(data, 1, hyperdis)
fish.result <- apply(data, 1, fishertest)
pval.f <- fish.result[1,]
est.f <- fish.result[2,]

data1 <- data
pvalue.mut  <- pvalue
pval.f.mut <- pval.f
est.f.mut <- est.f

df.mut <- fread("/abi/data/chen/project_august/mutation25_tial1/mut_pa.txt",data.table = F)

df.1 <- df.mut[,seq(1, ncol(df.mut), 4) ]
df.2 <- df.mut[,seq(2, ncol(df.mut), 4) ]
df.3 <- df.mut[,seq(3, ncol(df.mut), 4) ]
df.4 <- df.mut[,seq(4, ncol(df.mut), 4) ]

vt1 <- rowSums(df.1)
vt2 <- rowSums(df.2)
vt3 <- rowSums(df.3)
vt4 <- rowSums(df.4)

df.merge <- data.frame(cbind(vt1,vt2,vt3,vt4))
data <- df.merge
data2 <- data
pvalue <- apply(data, 1, hyperdis)
fish.result <- apply(data, 1, fishertest)
pval.f <- fish.result[1,]
est.f <- fish.result[2,]

pvalue.pa  <- pvalue
pval.f.pa <- pval.f
est.f.pa <- est.f





pvals <- data.frame(cbind(pvalue.mut,pvalue.pa))
pvalf <- data.frame(cbind(pval.f.mut,pval.f.pa))
estfs <- data.frame(cbind(est.f.mut,est.f.pa))


gc <- fread("/abi/data/chen/newhome/basic/coding.txt",data.table = F)
colnames(gc)[1] <- "chr"
gc <- gc[,c(1,4)]

pvals_tial1 <- cbind(gc,pvals)
pvals_f <- cbind(gc,pvalf)
ests_f <- cbind(gc,estfs)

all <- cbind(gc,data1,data2,pvals,pvalf,estfs)
names(all) <- c("chr","gene_name","nsample","mut_tial1","mut_gene","mut_both","n.2","mut_tial1.2","pa_gene","mut_pa","pv.mut","pv.pa",
                "pv.f.mut","pv.f.pa","est.f.mut","est.f.pa")
all$d1 <- all$mut_both-all$mut_gene
all$d2 <- all$mut_pa - all$pa_gene
