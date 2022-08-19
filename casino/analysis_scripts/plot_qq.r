library(ggplot2)
library(grid)
library(gridExtra)
library(directlabels)

options(echo=TRUE)
args<-commandArgs(trailingOnly=TRUE)

#coh<-"PRAD-US"
coh<-args[1] #the 1st argument is name of cohort

#analysis_path="/ibios/co02/chen/casino/results/analysis/"

real<-args[2] #the 2nd argument is file of scores
#real<-"/ibios/co02/chen/work/17_05_2016/enhancer/results/analysis/MB_snv_score.txt"
output_path<-args[3] #the 3rd argument is output PATH

#real<-paste(analysis_path,coh,"_score.bed",sep="")
a<-read.table(real,header=F,sep="\t",comment.char="")

#colnames(a)<-c("chr","start","end","gene","snv_score","cnv_score","length","mul_score")
colnames(a)<-c("chr","start","end","gene","snv_score")

pdf(paste(output_path,"/qqplot","_",coh,".pdf",sep=""))

for (i in 1:10){

#for each turn of randomization
#input file includes only 5 columns:chromosome,start,end,gene names and score.
#the file name of input should be ${coh}_score_snv${i}.bed
  
#rand <- paste("/ibios/co02/chen/random_10turns/results_promoter/analysis/",coh,"_score_snv",i,".bed",sep="")
rand <- paste("/ibios/co02/chen/ran_temp/MB_random/",coh,"_score_random_t",i,".txt",sep="")
b <- read.table(rand,header=F,sep="\t",comment.char="")
#b <- b[1,2,3,4,5]
colnames(b) <- c("chr","start","end","gene","snv_score")
df1 <- data.frame(a$chr,a$start,a$end,a$snv_score,a$gene)
df2 <- df1[ order(df1$a.snv_score), ]
colnames(df2)<-c("chr","start","end","real","gene")
df2$random<-sort(b$snv_score)
df2$var<-df2$real-df2$random
lims <- c(0,max(range(df2$real)[2],range(df2$random)[2]))
ns1 <- df2
#write.table(head(ns1[order(-ns1$var),],20),
#           file=paste("/home/hongc/test/qqplottest/qqplot_top_t",i,"_",coh,".txt",sep=""),
#           quote=F,row.names=F)
write.table(head(ns1[order(-ns1$var),],500),
            file=paste(output_path,"/qqplot",i,"_",coh,".txt",sep=""),
            quote=F,row.names=F,sep="\t")
#write.table(head(ns1[order(ns1$var),],20),
#           file=paste("/home/hongc/test/qqplottest/qqplot_bot_t",i,"_",coh,".txt",sep=""),
#           quote=F,row.names=F)
#print((ggplot(nsx,aes(x=random, y=real)) 
 #      + geom_point() + geom_abline(intercept = 0, slope = 1,color="red")
 #      + xlab('random mutation') + ylab('real data')+xlim(lims)+ylim(lims)+ggtitle(paste(coh,"_t",i,sep=""))
 #      + geom_dl(aes(label=ifelse(var>min(tail(sort(var),50)),as.character(gene),"")),list("smart.grid",cex=0.3))))
 
#geom_text(aes(label=Name),hjust=0, vjust=0)
}

dev.off()