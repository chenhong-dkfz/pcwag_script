#this code is used to calculate scores for SNV,INDEL,CNV and SV
#input file format:
#header should be with comment character "#"
#first four columns should be chromosome,start,end and gene name
#last line of CNV or SV input should be total length of patient, starting WITH "#"
#it is very important
#last line of SNV or INDEL input should be total length of patient, starting WITH "#"
#it is useless but necessary for file generated from MTATRIX_CONVERTER

options(echo=TRUE)
args<-commandArgs(trailingOnly = TRUE)

#################
#packages import#
#################

library(ggplot2)
library(grid)
library(gridExtra)
library(directlabels)

################
#Test arguments#
################
# arguments: input snv file, input cnv file, output score file
snv <- read.table(args[1],header=T,comment.char="")
cnv <- read.table(args[2],header=T,comment.char="")
#cnv <- read.table("/ibios/co02/chen/new_casino/matrix/promoter_matrix/SSS/PRAD-UK_cnv.bed",header=T,comment.char="")
#snv <- read.table("/ibios/co02/chen/new_casino/matrix/promoter_matrix/SSS/PRAD-UK_snv.bed",header=T,comment.char="")
ignored_chr_X <- 0

#snv scoring
snv <- snv[which(snv$X.chr!="X"),] #remove specialized chromosome information
snv <- snv[which(snv$X.chr!="Y"),]
#s0 <- snv[-nrow(snv),] #remove last line of cnv length information
s0 <- snv
ll <- dim(s0)[1]
s1 <- s0[-c(1,2,3,4)]
nopid<-dim(s1)[2]
#ss<-sapply(s1,function(x)-log((sum(x)-x+1)/(sum(x)+1),base=10))*(1/(nopid)) #redo it
ss<-sapply(s1,function(x)-log((x[ll]-x+1)/(x[ll]+1),base=10))*(1/(nopid))
ss2<-data.frame(ss)
ss2$sum <- rowSums(ss)
#normalize snv score, only for 1.0 version
#ss2$sum <-(ss2$sum-min(ss2$sum))/(max(ss2$sum)-min(ss2$sum))
snv_scored <- cbind(s0[c(1,2,3,4)],ss2)

#cnv-part#
##cnv scoring
#cnv <- cnv[which(cnv$X.chr!="X"),]
#cnv <- cnv[which(cnv$X.chr!="Y"),]
#genome_length = 3000000000
#c1 <- cnv[-c(1,2,3,4)] 
#ll <- dim(c1)[1] #line number of cnv total length
###cc<-sapply(c1,function(x)((-x[ll]+x-1)/(x[ll]+1))*log(((x+1)/genome_length),base=10)) # not any more
#cc<-sapply(c1,function(x)((-x-1)/(x[ll]+1))*log(((x+1)/genome_length),base=10))
#cc2 <- cc[-nrow(cc),] 
#cc2<-data.frame(cc2)
#cc2$sum <- rowSums(cc2)
##normalize cnv score, only for 1.0 version
##cc2$sum <-(cc2$sum-min(cc2$sum))/(max(cc2$sum)-min(cc2$sum))
#cnv_scored <- cbind(s0[c(1,2,3,4)],cc2)

##weight and output
#mix_score <- cnv_scored[c(1,2,3,4)]
#mix_score$snv_score <- snv_scored$sum
#mix_score$cnv_score <- cnv_scored$sum
#mix_score$length <- as.numeric(as.character(mix_score$end)) - as.numeric(as.character(mix_score$start))
#mix_score$snv_rs <- mix_score$snv_score * (1000/mix_score$length)
#mix_score$cnv_rs <- mix_score$cnv_score * (1000/mix_score$length)

#replace previous set#
mix_score <- snv_scored[c(1,2,3,4)]
mix_score$snv_score <- snv_scored$sum
mix_score <- mix_score[-nrow(mix_score),]
mix_score$length <- as.numeric(as.character(mix_score$end)) - as.numeric(as.character(mix_score$start))
mix_score$snv_rs <- mix_score$snv_score * (1000/mix_score$length)
#write.table(mix_score,args[3],sep="\t",quote=F,row.names=F,col.names=F)


coh<-args[5]


#random_matrix <- mix_score[c(1,2,3,4,8,9)]
random_matrix <- mix_score[c(1,2,3,4,6,7)]

for (i in 1:10){
  
  #for each turn of randomization
  #it is very similar to real input mutation matrix
  #but there is not last line for comment information
  #the file name of input should be ${coh}_t${i}_snv.bed
    
  rand <- paste(args[3],"/",coh,"/",coh,"_t",i,"_snv.bed",sep="")
  
  #rand <- paste("/ibios/co02/chen/new_casino/matrix/promoter_matrix/",coh,"/",coh,"_t",i,"_snvss.bed",sep="")
  b <- read.table(rand)
  
  #snv scoring
  b1 <- b[which(b$V1!="X"),] #remove specialized chromosome information
  b1 <- b1[which(b1$V1!="Y"),]
  
  names(b1)<-names(s0[dim(s0)[1],])
  b1 <- rbind(b1,s0[dim(s0)[1],])
  
  b2 <- b1[-c(1,2,3,4)]
  nopid<-dim(b2)[2]
  #bs<-sapply(b2,function(x)-log((sum(x)-x+1)/(sum(x)+1),base=10))*(1/(nopid))
  bs<-sapply(b2,function(x)-log((x[ll]-x+1)/(x[ll]+1),base=10))*(1/(nopid))
  bs2<-data.frame(bs)
  bs2$sum <- rowSums(bs)
  #normalize snv score, only for 1.0 version
  #ss2$sum <-(ss2$sum-min(ss2$sum))/(max(ss2$sum)-min(ss2$sum))
  bs_scored <- cbind(b1[c(1,2,3,4)],bs2)
  bs_scored <- bs_scored[-nrow(bs_scored),]
  
  bs_scored$length <- as.numeric(as.character(bs_scored[,3])) - as.numeric(as.character(bs_scored[,2]))
  bs_scored$snv_rs <- bs_scored$sum * (1000/bs_scored$length)
  
  df1 <- data.frame(mix_score$X.chromosome,mix_score$start,mix_score$end,mix_score$gene,mix_score$snv_rs)
  df1$cd <- 1:dim(df1)[1]
  df2 <- df1[ order(df1$mix_score.snv_rs), ]
  colnames(df2)<-c("chr","start","end","gene","real","cd")
  df2$random<-sort(bs_scored$snv_rs)
  df2$var<-df2$real-df2$random
  df3 <- df2[ order(df2$cd), ]
  random_matrix[,paste("var",i,sep="")]<-as.numeric(as.character(df3$var))
  
  
}

for (i in 1:10){
  var_threshold <- sort(random_matrix[,paste("var",i,sep="")],decreasing=T)[20]
  random_matrix[,paste("hl",i,sep="")] <- ifelse(random_matrix[,paste("var",i,sep="")]>=var_threshold,1,0)
  }

random_matrix2 <- random_matrix
#qqdraw <- random_matrix2[,c(4,5,7:16)]
random_matrix$hls <- rowSums(random_matrix[,17:26])
#random_matrix <- random_matrix[,c(1,2,3,4,6,27)]

write.table(random_matrix,file=args[4],row.names=F,quote=F,col.names=T)
#for (i in 1:10){
#  qqtmps <- qqdraw[ order(-qqdraw[,paste("var",i,sep="")]), ]
#  qqtmp <- head(qqtmps,100000)
#  lims <- c(0,max(range(qqtmp$snv_rs)[2],range(qqdraw[,paste("var",i,sep="")])[2]))
#  p <- (ggplot(qqtmp,aes(x=qqtmp[,paste("var",i,sep="")],y=snv_rs))
#        +geom_point()+geom_abline(intercept=0,slope=1,color="red")
#        +xlab("random_mutation")+ylab("real data")+xlim(lims)+ylim(lims)+ggtitle(paste("turn",i,sep="_")))
#  p
#}