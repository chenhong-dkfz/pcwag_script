options(echo=TRUE)
args<-commandArgs(trailingOnly = TRUE)

#cnv <- read.table("/ibios/co02/chen/enhancer/results/analysis/EOPC-DE_cnv.bed",header=T,comment.char="")
#snv <- read.table("/ibios/co02/chen/enhancer/results/analysis/EOPC-DE_snv.bed",header=T,comment.char="")

#######################
#here is only for test#
#######################
#snv <- read.table("/home/hongc/examples/test_snv.bed",header=T,comment.char="")




#cnv <- read.table(args[1],header=T,comment.char="")
snv <- read.table(args[1],header=T,comment.char="",fill=T)

#snv scoring
colnames(snv)[1] <- "chromosome"
#snv <- snv[which(snv$chromosome!="X"),] 
#s0 <- snv[-nrow(snv),]
s0 <- snv
s1 <- s0[-c(1,2,3,4)]
nopid<-dim(s1)[2]
ll <- dim(s1)[1]
#ss<-sapply(s1,function(x)-log((x[ll]-x+1)/(x[ll]+1),base=10))*(1/(nopid))
ss<-sapply(s1,function(x)-log((sum(x)-x+1)/(sum(x)+1),base=10))*(1/(nopid))
ss2<-data.frame(ss)
ss2$sum <- rowSums(ss)
snv_scored <- cbind(s0[c(1,2,3,4)],ss2$sum)
colnames(snv_scored)[5] <- "snv_score"
#snv_scored<-snv_scored[-dim(snv_scored)[1],]
#args[3]
write.table(snv_scored,args[2],sep="\t",quote=F,row.names=F,col.names=F)
