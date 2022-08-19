options(echo=TRUE)
args<-commandArgs(trailingOnly = TRUE)


#######################
#here is only for test#
#######################
#snv <- read.table("/home/hongc/examples/test_snv.bed",header=T,comment.char="")



#read files
cnv <- read.table(args[1],header=T,comment.char="")
#snv <- read.table(args[2],header=T,comment.char="")

#snv scoring
#snv <- snv[which(snv$X.chr!="X"),]
#s0 <- snv[-nrow(snv),]
#s1 <- s0[-c(1,2,3,4)]
#nopid<-dim(s1)[2]
#ss<-sapply(s1,function(x)-log((sum(x)-x+1)/(sum(x)+1),base=10))*(1/(nopid))
#ss2<-data.frame(ss)
#ss2$sum <- rowSums(ss)
##normalize snv score, only for 1.0 version
##ss2$sum <-(ss2$sum-min(ss2$sum))/(max(ss2$sum)-min(ss2$sum))
#snv_scored <- cbind(s0[c(1,2,3,4)],ss2)

#cnv scoring
cnv <- cnv[which(cnv$X.chr!="X"),]
s0 <- cnv[-nrow(cnv),]
genome_length = 3000000000
c1 <- cnv[-c(1,2,3,4)]
ll <- dim(c1)[1] #cnv total length
#cc<-sapply(c1,function(x)((-x[ll]+x-1)/(x[ll]+1))*log(((x+1)/genome_length),base=10))
cc<-sapply(c1,function(x)((-x-1)/(x[ll]+1))*log(((x+1)/genome_length),base=10))
cc2 <- cc[-nrow(cc),]
cc2<-data.frame(cc2)
cc2$sum <- rowSums(cc2)
#normalize cnv score, only for 1.0 version
#cc2$sum <-(cc2$sum-min(cc2$sum))/(max(cc2$sum)-min(cc2$sum))
cnv_scored <- cbind(s0[c(1,2,3,4)],cc2)

#weight and output
#w_snv = 100
#w_cnv = 1
mix_score <- cnv_scored[c(1,2,3,4)]
#mix_score$snv_score <- snv_scored$sum
mix_score$cnv_score <- cnv_scored$sum
colnames(mix_score)<-c("chr","start","end","gene","cnv_score")
mix_score$length <- as.numeric(as.character(mix_score$end)) - as.numeric(as.character(mix_score$start))+1
#final score calculated by sum, only for 1.0 version
##mix_score$score <- (w_snv*snv_scored$sum + w_cnv*cnv_scored$sum)*(1000/(as.numeric(as.character(mix_score$end)) - as.numeric(as.character(mix_score$start))+1))
##mix_score$mul_score <- (w_snv*snv_scored$sum * w_cnv*cnv_scored$sum)
mix_score$fs <- as.numeric(as.character(mix_score$cnv_score))/as.numeric(as.character(mix_score$length))
nd <- mix_score[order(-mix_score$fs),]
write.table(nd,args[2],sep="\t",quote=F,row.names=F,col.names=F)

