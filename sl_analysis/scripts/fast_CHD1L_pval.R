
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
  res <- fisher.test(mtx,alternative = "less")
  p.value <- res$p.value
  odds <- res$estimate[[1]]
  return(c(p.value,odds))
}

mutss <- c("mut_mut","mut_nam","mut_pam","mut_pa","mut_ra",
           "nam_mut","nam_nam","nam_pam","nam_pa","nam_ra",
           "pam_mut","pam_nam","pam_pam","pam_pa","pam_ra",
           "ra_mut","ra_nam","ra_pam","ra_pa","ra_ra",
           "pa_mut","pa_nam","pa_pam","pa_pa","pa_ra")

gene.name <- "CHD1L"

# file input : n_total, n_mut1, n_mut2, n_intersect
# comet input, c(n_total,n_mut1, n_mut2, n_intersect)
# phyper input, n_intersect, n_mut1, n_total-n_mut1, n_mut2

pvalue <- list()
types <- list()
pval.f <- list()
est.f <- list()

for (i in 1:length(mutss)){
  muts <- mutss[i]
  print(muts)
  data <- fread(paste0("/abi/data/chen/project_august/mutation25/cohortwise/",gene.name,"/",muts,".txt"),
                data.table=F)
  #data <- data[1:100,]
  pvalue[[i]] <- apply(data, 1, hyperdis)
  fish.result <- apply(data, 1, fishertest)
  pval.f[[i]] <- fish.result[1,]
  est.f[[i]] <- fish.result[2,]
}

pvals <- data.frame(matrix(unlist(pvalue),ncol = 25))
pvalf <- data.frame(matrix(unlist(pval.f),ncol = 25))
estfs <- data.frame(matrix(unlist(est.f),ncol = 25))

colnames(pvals) <- mutss
colnames(pvalf) <- mutss
colnames(estfs) <- mutss

gc <- fread("/home/hongc/PhD/basic/coding.txt",data.table = F)
colnames(gc)[1] <- "chr"
gc <- gc[,c(1,4)]

pvals_chd1l <- cbind(gc,pvals)
pvals_f <- cbind(gc,pvalf)
ests_f <- cbind(gc,estfs)



write.table(pvals_chd1l,file="/abi/data/chen/project_august/mutation25/cohortwise/full2/alc1_hyper.txt",
            col.names = T, row.names = F, quote = F, sep = "\t")
write.table(pvals_f,file="/abi/data/chen/project_august/mutation25/cohortwise/full2/alc1_fisher.txt",
            col.names = T, row.names = F, quote = F, sep = "\t")
write.table(ests_f,file="/abi/data/chen/project_august/mutation25/cohortwise/full2/alc1_oddsratio.txt",
            col.names = T, row.names = F, quote = F, sep = "\t")
