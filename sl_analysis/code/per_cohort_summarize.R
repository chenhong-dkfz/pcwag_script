library(data.table)

cohorts = c("BLCA-US","CLLE-ES","GBM-US","LGG-US","MALY-DE","PAEN-AU","READ-US",
           "BOCA-UK","CMDI-UK","HNSC-US","LICA-FR","MELA-AU","PBCA-DE","RECA-EU",
           "BRCA-EU","COAD-US","KICH-US","LIHC-US","ORCA-IN","SARC-US",
           "BRCA-UK","DLBC-US","KIRC-US","LINC-JP","OV-AU","SKCM-US",
           "BRCA-US","EOPC-DE","KIRP-US","LIRI-JP","OV-US","PRAD-CA","STAD-US",
           "BTCA-SG","ESAD-UK","LAML-KR","LUAD-US","PACA-AU","PRAD-UK","THCA-US",
           "CESC-US","GACA-CN","LAML-US","LUSC-US","PACA-CA","PRAD-US","UCEC-US")
m.cohort <- cohorts[c(28,32,39,46)] # prostate
m.cohort <- c("BRCA-EU","BRCA-UK","BRCA-US") # breast cancer
m.cohort <- c("SKCM-US","MELA-AU")
m.cohort <- c("KICH-US")
cohort.hist <- fread("/home/hongc/PhD/basic/cohort_histology.tsv",data.table = F,header = F)




muts <- c("mut","nam","pam","pa","ra")
mutss <- NULL
i <- 1
for (mut in muts){
  for (mut2 in muts){
    st <- paste0(mut,"_",mut2)
  mutss[i] <- st
  i <- i + 1
  }
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



gene_list <- fread("/home/hongc/PhD/basic/coding.txt",data.table = F)
gene_list$id <- rownames(gene_list)
t1 <- fread("/home/hongc/test/t1.txt",data.table = F)
t2 <- fread("/home/hongc/Downloads/gene_list_PARPInhSens.tsv",data.table = F)
t1$PARPInhSens <- FALSE
t1[t1$gene_name %in% t2$`Associated Gene Name`,4] <- TRUE
genes <- cbind(t1,gene_list)
gene_list <- genes[,c(1,6,7,2,3,4)]

# per gene #


gene <- "SMAD4"
gene_id <- gene_list[gene_list$gene_name==gene,5]

meta <- NULL

for(i in 1:length(m.cohort)){
  
  coh.4v <- NULL
  for (mut in mutss){
  print(mut)
  fin <- fread(paste0("/abi/data/chen/project_august/mutation25/per_cohort/",m.cohort[i],"/PTEN2/",mut,".txt"),data.table = F)
  coh.4v <- rbind(coh.4v,fin[gene_id,])
  }
  rownames(coh.4v)<- mutss
  meta[[i]] <- coh.4v
}

meta.cohort <- Reduce("+",meta)
fish.result <- apply(meta.cohort, 1, fishertest)
meta.cohort$p <- fish.result[1,]
meta.cohort$odds <- fish.result[2,]

colnames(meta.cohort) <- c("sample_size","PTEN_mut",paste0(gene,"_mut"),"overlaps","p","OR")
meta.cohort$mutation_type <- rownames(meta.cohort)
meta.cohort <- meta.cohort[,c(7,1,2,3,4,5,6)]
write.table(meta.cohort,file=paste0("/home/hongc/test/pten_",gene,"_prostate_cancer.txt"),sep = "\t",quote = F,row.names = F)


#save(meta.pv,file="/abi/data/chen/project_august/mutation25/per_cohort/p_value_PTEN.rdata")
#save(meta.odd,file="/abi/data/chen/project_august/mutation25/per_cohort/odds_ratio_PTEN.rdata")

# for each gene #
cohort.hist$V2 <- gsub("/", "_", cohort.hist$V2)
cohort.hist <- cohort.hist[cohort.hist$V1!="PAEN-IT",]
histologies <- names(table(cohort.hist$V2))



for(hist in histologies){
  projects <- c(cohort.hist[cohort.hist$V2==hist,1])

for (mut in mutss){  
  
meta <- NULL
m.cohort <- projects
print(hist)
for(i in 1:length(m.cohort)){
  
  coh.4v <- NULL
  #mut <- "pa_pa"
  print(mut)
  fin <- fread(paste0("/abi/data/chen/project_august/mutation25/per_cohort/",m.cohort[i],"/CHD1L/",mut,".txt"),data.table = F)
  #fin <- fread(paste0("/abi/data/chen/project_august/mutation25/per_cohort/",m.cohort[i],"/PTEN2/",mut,".txt"),data.table = F)
  if(is.null(meta)==TRUE)
  {meta <-fin}else{
    meta <- meta+fin
  }

}


fish.result <- apply(meta, 1, fishertest)
meta$p <- fish.result[1,]
meta$odds <- fish.result[2,]
meta <- cbind(gene_list,meta)
colnames(meta) <- c("chr","start","end","gene","cgc","PARPInhSens","sample_size","ALC1_mut","gene_mut","overlaps","p","OR")
write.table(meta,file=paste0("/home/hongc/test/result_subcohort3/",hist,"_",mut,".txt"), sep = "\t",quote = F,row.names = F)
}
}
