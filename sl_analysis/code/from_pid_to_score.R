library(data.table)
library(ggplot2)


data <- fread("/abi/data/chen/project_august/mutation2.csv",data.table=F)


data.new <- data[data$gene_name=="PTEN"|data$gene_name=="CHD1L",]
rownames(data.new) <- data.new$gene_name
data.new <- data.new[,-1]
#data.new <- data.new[ , -which(names(data.new) %in% c("V444","V865"))]

data.new <- replace(data.new,data.new==5,"real biallelic")
data.new <- replace(data.new,data.new==4,"possible biallelic")
data.new <- replace(data.new,data.new==3,"monoallelic")
data.new <- replace(data.new,data.new==2,"AMP>8")
data.new <- replace(data.new,data.new==1,"AMP_4-8")
data.new <- replace(data.new,data.new==0,"")

mat <- data.new

cohort.data <- fread("/abi/data/chen/project_august/pcawg_ta_project.txt",data.table = F)
t_mat <- data.frame(t(mat))
t_mat$pid <- rownames(t_mat)
t_mat_s <- merge(t_mat,cohort.data,by.x="pid",by.y="tumor_aliquot_id_wgs")

#t_mat_s <- t_mat_s[t_mat_s$PTEN!="",]
p <- ggplot(t_mat_s,aes(x=PTEN,y=1,fill=PTEN))+geom_bar(stat="identity" )+facet_wrap(~project_code)
p <- p + theme(axis.title.y=element_blank(),axis.text.x = element_text(angle = -90))
# pdf("/home/hongc/test/result_100918/PTEN.pdf",width = 20)
# print(p)
# dev.off()


#
t_mat_s$project_size <- apply(t_mat_s, 1, function(x) table(t_mat_s$project_code)[x[4]])
t_mat_s$project_mono <- apply(t_mat_s, 1, function(x) table(t_mat_s$project_code)[x[4]])

chd1l.dis <- data.frame(table(t_mat_s$CHD1L,t_mat_s$project_code))
names(chd1l.dis) <- c("mutations","cohort","Freq")
chd1l.dis[chd1l.dis$mutations=="",]$mutations <- NA
p1 <- ggplot(chd1l.dis,aes(y=Freq,x=cohort,fill=mutations))+coord_flip() +
  geom_bar(stat = "identity",position = position_stack(reverse = TRUE))+
  ylab("Frequency")+xlab("cohort")+ggtitle("CHD1L mutation distributions")
# pdf("/home/hongc/test/result_100918/CHD1L_distribution.pdf",width = 20)
# print(p1)
# dev.off()

pten.dis <- data.frame(table(t_mat_s$PTEN,t_mat_s$project_code))
names(pten.dis) <- c("mutations","cohort","Freq")
pten.dis[pten.dis$mutations=="",]$mutations <- NA
p2 <- ggplot(pten.dis,aes(y=Freq,x=cohort,fill=mutations))+coord_flip() +
  geom_bar(stat = "identity",position = position_stack(reverse = TRUE))+
  ylab("Frequency")+xlab("cohort")+ggtitle("PTEN mutation distributions")
# pdf("/home/hongc/test/result_100918/PTEN_distribution.pdf",width = 20)
# print(p2)
# dev.off()

chd1l.dis <- data.frame(table(t_mat_s$CHD1L,t_mat_s$project_code))
names(chd1l.dis) <- c("mutations","cohort","Freq")
mylist <- split(chd1l.dis, chd1l.dis$mutations)
mylist_neg <- mylist[[1]]+mylist[[4]]
mylist_neg$cohort <- mylist[[1]]$cohort
mylist_pos <- mylist[[2]]+mylist[[3]]
mylist_pos$cohort <- mylist[[2]]$cohort
mylist_all <- cbind(mylist_pos,mylist_neg)
names(mylist_all) <- c("mut1","cohort1","freq1","mut2","cohort2","freq2")
mylist_all$freq_all <- mylist_all$freq1 + mylist_all$freq2
mylist_all$perc <- mylist_all$freq1/mylist_all$freq_all
chd1l.prop <- mylist_all
#save(mylist_all,file="/home/hongc/test/result_240918/chd1l.rdata")

pten.dis <- data.frame(table(t_mat_s$PTEN,t_mat_s$project_code))
names(pten.dis) <- c("mutations","cohort","Freq")
mylist <- split(pten.dis, pten.dis$mutations)
mylist_neg <- mylist[[1]]+mylist[[2]]
mylist_neg$cohort <- mylist[[1]]$cohort
mylist_pos <- mylist[[3]]+mylist[[4]]+mylist[[5]]
mylist_pos$cohort <- mylist[[3]]$cohort
mylist_all <- cbind(mylist_pos,mylist_neg)
names(mylist_all) <- c("mut1","cohort1","freq1","mut2","cohort2","freq2")
mylist_all$freq_all <- mylist_all$freq1 + mylist_all$freq2
mylist_all$perc <- mylist_all$freq1/mylist_all$freq_all
pten.prop <- mylist_all
#save(pten.prop,file="/home/hongc/test/result_240918/pten.rdata")


pten.top <- pten.prop[pten.prop$perc>0.15,c(2,3,6,8)]
chd1l.top <- chd1l.prop[chd1l.prop$perc>0.15,c(2,3,6,8)]

cohort.hist <- fread("/home/hongc/PhD/basic/cohort_histology.tsv",data.table = F,header = F)
names(cohort.hist) <- c("cohort","hist")
pten.top2 <- merge(pten.top,cohort.hist,by.x="cohort1",by.y="cohort")
chd1l.top2 <- merge(chd1l.top,cohort.hist,by.x="cohort1",by.y="cohort")
cohort.data2 <- merge(cohort.data,cohort.hist,by.x="project_code",by.y = "cohort")

# test version #

pten.top <- pten.prop[pten.prop$perc>-1,c(2,3,6,8)]
chd1l.top <- chd1l.prop[chd1l.prop$perc>-1,c(2,3,6,8)]
cohort.hist <- fread("/home/hongc/PhD/basic/cohort_histology.tsv",data.table = F,header = F)
names(cohort.hist) <- c("cohort","hist")
pten.top2 <- merge(pten.top,cohort.hist,by.x="cohort1",by.y="cohort")
chd1l.top2 <- merge(chd1l.top,cohort.hist,by.x="cohort1",by.y="cohort")
cohort.data2 <- merge(cohort.data,cohort.hist,by.x="project_code",by.y = "cohort")

write.table(pten.top2,file = "/home/hongc/test/result_240918/pten_cohort_dist.txt",
             quote = F, row.names = F, sep = "\t")

write.table(chd1l.top2,file = "/home/hongc/test/result_240918/chd1l_cohort_dist.txt",
            quote = F, row.names = F, sep = "\t")

# old version #

#############
# pten part #
#############

the.pid <-cohort.data[which(cohort.data$project_code %in% pten.top2$cohort1),] 
data.pten <- data[,colnames(data) %in% c(the.pid$tumor_aliquot_id_wgs,"gene_name")]
write.table(data.pten,file="/abi/data/chen/project_august/mutation_pten.csv",col.names = T,
            row.names = F, quote = F, sep = "\t")
##############
# chd1l part #
##############
the.pid <-cohort.data[which(cohort.data$project_code %in% chd1l.top2$cohort1),] 
data.chd1l<- data[,colnames(data) %in% c(the.pid$tumor_aliquot_id_wgs,"gene_name")]
write.table(data.pten,file="/abi/data/chen/project_august/mutation_chd1l.csv",col.names = T,
            row.names = F, quote = F, sep = "\t")

# old version ends #

# new version #

#############
# pten part #
#############

the.pid <-cohort.data2[which(cohort.data2$hist %in% pten.top2$hist),] 
data.pten <- data[,colnames(data) %in% c(the.pid$tumor_aliquot_id_wgs,"gene_name")]
write.table(data.pten,file="/abi/data/chen/project_august/mutation_pten.csv",col.names = T,
            row.names = F, quote = F, sep = "\t")
##############
# chd1l part #
##############
the.pid <-cohort.data2[which(cohort.data2$hist %in% chd1l.top2$hist),] 
data.chd1l<- data[,colnames(data) %in% c(the.pid$tumor_aliquot_id_wgs,"gene_name")]
write.table(data.chd1l,file="/abi/data/chen/project_august/mutation_chd1l.csv",col.names = T,
            row.names = F, quote = F, sep = "\t")

library(data.table)
library(ComplexHeatmap)
library(dplyr)
library(easyGgplot2)
library(grid)
library(gridExtra)

library(cometExactTest)
library(data.table)

#args = commandArgs(trailingOnly=TRUE)



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




#gene.name <- args[1]
#muts <- args[2]



#setwd("/abi/data")

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
  data <- data[1:100,]
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
pvals_chd1l <- pvals_chd1l[,c(1,2,8:17)]


gene.name <- "PTEN"

# file input : n_total, n_mut1, n_mut2, n_intersect
# comet input, c(n_total,n_mut1, n_mut2, n_intersect)
# phyper input, n_intersect, n_mut1, n_total-n_mut1, n_mut2

pvalue <- list()
types <- list()

for (i in 1:length(mutss)){
  muts <- mutss[i]
  print(muts)
  data <- fread(paste0("/abi/data/chen/project_august/mutation25/cohortwise/",gene.name,"/",muts,".txt"),
                data.table=F)
  
  pvalue[[i]] <- apply(data, 1, hyperdis)
  
}

pvals <- data.frame(matrix(unlist(pvalue),ncol = 25))
colnames(pvals) <- mutss
#gc <- fread("/home/hongc/PhD/basic/coding.txt",data.table = F)
#colnames(gc)[1] <- "chr"
#gc <- gc[,c(1,4)]
pvals_pten <- cbind(gc,pvals)
pvals_pten <- pvals_pten[,c(1,2,3:7,18:27)]

write.table(pvals_chd1l,file="/abi/data/chen/project_august/mutation25/cohortwise/pvals/chd1l.txt",
            col.names = T, row.names = F, quote = F, sep = "\t")
write.table(pvals_pten,file="/abi/data/chen/project_august/mutation25/cohortwise/pvals/pten.txt",
            col.names = T, row.names = F, quote = F, sep = "\t")



#qvals_chd1l_body <- sapply(pvals_chd1l[,3:12],p.adjust,method="fdr")
#qvals_chd1l <- cbind(pvals_chd1l[,1:2],qvals_chd1l_body)

pvals_chd1l <- fread("/abi/data/chen/project_august/mutation25/cohortwise/pvals/chd1l.txt",data.table=F)
pvals_pten <- fread("/abi/data/chen/project_august/mutation25/cohortwise/pvals/pten.txt",data.table=F)
cgc.table <- fread("/home/hongc/PhD/basic/cgc.txt",data.table = F)
cgc.table <- cgc.table[cgc.table$`Role in Cancer`!="",]
chd1l.selected <- merge(pvals_chd1l,cgc.table,by.x="gene_name",by.y="Gene Symbol")
pten.selected <- merge(pvals_pten,cgc.table,by.x="gene_name",by.y="Gene Symbol")
write.table(pten.selected,file="/home/hongc/test/result_240918/pten_cgc.txt",sep = "\t",row.names = F,quote = F)
write.table(chd1l.selected,file="/home/hongc/test/result_240918/chd1l_cgc.txt",sep = "\t",row.names = F,quote = F)
