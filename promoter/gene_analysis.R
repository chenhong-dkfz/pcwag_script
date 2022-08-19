suppressMessages(library(data.table))
#suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(reshape2))
suppressMessages(library(grid))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(gtools))
suppressMessages(library(entropy))
suppressMessages(library(lattice))

# promoter analysis

# gene names

args = commandArgs(trailingOnly=TRUE)
gene.name <- args[1]
output.dir <- args[2]
#gene.name <- "PTEN"


# read files

#snv.table <- fread("/home/hongc/PhD/TumorPrint/single_gene/test/snv.txt",data.table = F)
#cnv.table <- fread("/home/hongc/PhD/TumorPrint/single_gene/test/cnv.txt",data.table = F)
#exp.table <- fread("/home/hongc/PhD/TumorPrint/single_gene/test/exp.txt",data.table = F)

snv.table <- fread("/ibios/co02/chen/pcawg_2018/result/matrix/promoter/report/pm_top.csv",data.table = F)
exp.table <- fread("/ibios/co02/chen/pcawg_2018/result/matrix/promoter/report/exp_file_ta.csv",data.table = F)

#blacklist <- c("AQP1","CRYBG3","IER3IP1","IGHV3OR16-13","MRPS17","SIGLEC5","SNORA40","SNURF")

#'%!in%' <- function(x,y)!('%in%'(x,y))
#snv.table <- snv.table[snv.table$gene %!in% blacklist,]

#colnames(exp.table)[1] <- "gene"
#exp.table <- exp.table[exp.table$gene %!in% blacklist,]

rownames(snv.table) <- snv.table[, 4] 
snv.table <- snv.table[, -c(1,2,3,4)]  
#rownames(cnv.table) <- cnv.table[, 1] 
#cnv.table <- cnv.table[, -1]  
rownames(exp.table) <- exp.table[, 1] 
exp.table <- exp.table[, -1]  

####################
# definition tables#
####################

pid.cohort <- fread("/ibios/co02/chen/pcawg_2018/pid_list/pid.tsv",data.table = F,header = F)
rownames(pid.cohort) <- pid.cohort[, 1] 
#pid.cohort <- data.frame(c("P1","P2","P3","P4"),c("C1","C1","C2","C2"))
colnames(pid.cohort) <- c("pid","cohort")
#rownames(pid.cohort) <- pid.cohort[,1]

# combine tables
cohort.meta <- fread("/home/hongc/PhD/basic/cohort_histology.tsv",data.table = F,header = F)
colnames(cohort.meta) <- c("cohort","meta")
pid.meta <- merge(pid.cohort,cohort.meta,by="cohort")

genes <- rownames(snv.table)

p.table <- data.frame(rep(0,length(genes)),rep(0,length(genes)),rep(0,length(genes)),rep(0,length(genes)))
colnames(p.table) <- c("gene","pvalue","exp_wt","exp_mt")
p.table$gene <- genes

snv.t <- data.frame(t(snv.table),check.names = F)
exp.t <- data.frame(t(exp.table),check.names = F)
snv.t$pid <- rownames(snv.t)
exp.t$pid <- rownames(exp.t)

for(gene in genes){
print(gene)
genei <- which(colnames(snv.t) %in% gene)  
genej <- which(colnames(exp.t) %in% gene)  
snv.t2 <- snv.t[,c(genei,ncol(snv.t))]
exp.t2 <- exp.t[,c(genej,ncol(exp.t))]
gene.matrix <- merge(snv.t2,exp.t2,by="pid")
gene.matrix <- na.omit(gene.matrix)

if (nrow(gene.matrix)== 0){
  p.table[genei,2] <- "NA"
  p.table[genei,3] <- "NA"
  p.table[genei,4] <- "NA"
} else if(length(table(gene.matrix[,2])) == 1){
  p.table[genei,2] <- "NA"
  p.table[genei,3] <- "NA"
  p.table[genei,4] <- "NA"
} else {
colnames(gene.matrix) <- c("pid","snv","exp")
gene.matrix$snv <- as.factor(gene.matrix$snv)

g0 <- gene.matrix[gene.matrix$snv==0,]
g1 <- gene.matrix[gene.matrix$snv==1,]
exp0 <- median(g0$exp)
exp1 <- median(g1$exp)

p0 <- wilcox.test(gene.matrix$exp~gene.matrix$snv,alternative="less")$p.value
#ggplot(gene.matrix,aes(y=exp,x=snv))+geom_boxplot()

p.table[genei,2] <- p0
p.table[genei,3] <- exp0
p.table[genei,4] <- exp1
}
}

p.table <- p.table[order(p.table$pvalue),]
ps <- p.table[p.table$pvalue<0.05,]
ps$exp_wt <- as.numeric(as.character(ps$exp_wt))
ps$exp_mt <- as.numeric(as.character(ps$exp_mt))
ps$diff <- ps$exp_mt - ps$exp_wt

write.table(p.table,file="/ibios/co02/chen/pcawg_2018/result/matrix/promoter/report/pv10k.txt",
            quote = F,row.names = F,col.names = T,sep="\t")

sgs <- ps$gene
#sg <- "PLA2G4A"

pdf("/home/hongc/test/propm.pdf",width = 15)

for (sg in sgs){
  
  
sg <- "C1orf112"  
genei <- which(colnames(snv.t) %in% sg)  
genej <- which(colnames(exp.t) %in% sg)  
snv.t2 <- snv.t[,c(genei,ncol(snv.t))]
exp.t2 <- exp.t[,c(genej,ncol(exp.t))]
gene.matrix <- merge(snv.t2,exp.t2,by="pid")
gene.matrix <- na.omit(gene.matrix)
colnames(gene.matrix) <- c("pid","snv","exp")
gene.matrix$snv <- as.factor(gene.matrix$snv)
gene.matrix2 <- merge(gene.matrix,pid.meta,by="pid")

give.n <- function(x){
  return(c(y = 7, label = length(x)))
}

frequencyAnnotation <- function(x) {
  c(y = (quantile(x, .75, names = F) + median(x))/2 + 2 , label=length(x))
}

pi0 <- ggplot(gene.matrix2,aes(y=exp,x=meta,fill=snv))+geom_boxplot()+
  stat_summary(fun.data = frequencyAnnotation, geom='text', 
               position = position_dodge(width = 0.75)) +
  ggtitle(paste0(sg," point mutations in promoters"))
plot(pi0)



}

dev.off()


pdf(paste0("/home/hongc/PhD/promoters/result/",sg,".pdf"),width = 15)
plot(pi0)
dev.off()
