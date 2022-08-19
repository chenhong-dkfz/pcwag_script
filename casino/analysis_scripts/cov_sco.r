scores <- read.table("/ibios/co02/chen/casino/gencode/results/analysis/MB_score.bed",header=F)
colnames(scores) <- c("chromosome","start","end","gene","score")
rt <- read.table("/ibios/co02/chen/casino/gencode/results/rep2_uniq.bed",header=F)
colnames(rt) <- c("chromosome","start","end","gene","rt")
hic <- read.table("/ibios/co02/chen/casino/gencode/results/hic2_uniq.bed",header=F)
colnames(hic) <- c("chromosome","start","end","gene","hic")

rt2 <- rt[which(rt$chromosome!="X"),]
rt2 <- rt2[-dim(rt2)[1],]
hic2 <- hic[which(hic$chromosome!="X"),]
hic2 <- hic2[-dim(hic2)[1],]

scores2 <- scores
scores2$rt <- rt2$rt
scores2$hic <- hic2$hic

rt0 <- scores2[which(scores2$rt!=0),]
hic0 <- scores2[which(scores2$hic!=0),]


rt0 <- rt0[which(rt0$score!=0),]
p <- ggplot(rt0, aes(rt,score))
p1 <- p + geom_point(aes(colour = factor(chromosome)), alpha = 0.6, position = 'jitter')
print(p1)

raw <- read.table("/ibios/co02/chen/casino/gencode/results/analysis/MB_snv.bed",header=T)
s1 <- raw[-c(1,2,3,4)]
s1$sum <- rowSums(s1)
