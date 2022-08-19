library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)

#import data

setwd("/abi/data/chen/newhome/promoter_focallity/data/")
scores <- fread("intermediate/int2.tsv",data.table = F)
counts1 <- fread("all_counts_1e6.txt",data.table = F)
counts2 <- fread("all_counts_1e7.txt",data.table = F)

#gene.list <- fread("/ibios/co02/chen/pcawg_2018/result/matrix/coding/LGG-US2.bed",data.table = F)
gene.list <- fread("/abi/data/chen/co02/pcawg_2018/result/matrix/coding/LGG-US2.bed",data.table = F)
gene.list <- gene.list[,c(1,2,3,4)]
names(gene.list) <- c("chr","start","end","gene")

names(scores) <- c("fscore_1e6","fscore_1e7","nb_1e6","nb_1e7","local_1e6","local_1e7","peak_1e6","peak_1e7")
scores$counts_1e6 <- counts1$V1
scores$counts_1e7 <- counts2$V1

df <- cbind(gene.list,scores)

cgc <- fread("/abi/data/chen/newhome/basic/cgc.txt",data.table = F)
df$cgc <- "none"
df[df$gene%in%cgc$`Gene Symbol`,"cgc"] <- "cgc"
df$tsg <- "none"
df[df$gene%in%cgc[cgc$`Role in Cancer`%like%"TSG","Gene Symbol"],"tsg"] <- "tsg"

df$fs <- "none"
fs.list <- fread("/abi/data/chen/newhome/fragile_site_bed/gene.txt",data.table = F)
names(fs.list) <- c("chr","start","end","gene","score","strand")
df[df$gene%in%fs.list$`gene`,"fs"] <- "fragile"


# focallity scores
compare_means(fscore_1e6 ~ tsg, data = df)
p <- ggplot(data = df, aes(y=fscore_1e6,x=tsg,fill=tsg))+geom_boxplot(outlier.shape = NA)
ylim1 = boxplot.stats(df$fscore_1e6)$stats[c(1, 5)]
p_fs_1e6 = p + coord_cartesian(ylim = ylim1*5) + stat_compare_means(label.x = 0.6, label.y = 6)+ scale_fill_manual(values=c("deepskyblue4", "firebrick2"))
p_fs_1e6

compare_means(fscore_1e7 ~ tsg, data = df)
p <- ggplot(data = df, aes(y=fscore_1e7,x=tsg,fill=tsg))+geom_boxplot(outlier.shape = NA)
ylim1 = boxplot.stats(df$fscore_1e7)$stats[c(1, 5)]
p_fs_1e7 = p + coord_cartesian(ylim = ylim1*5) + stat_compare_means(label.x = 0.6, label.y = 50)+ scale_fill_manual(values=c("deepskyblue4", "firebrick2"))
p_fs_1e7

# nb ranks
compare_means(nb_1e6 ~ tsg, data = df)
p <- ggplot(data = df, aes(y=nb_1e6,x=tsg,fill=tsg))+geom_boxplot(outlier.shape = NA)
ylim1 = boxplot.stats(df$nb_1e6)$stats[c(1, 5)]
p_nb_1e6 = p + coord_cartesian(ylim = ylim1*1.5) + stat_compare_means(label.x = 0.6, label.y = 12)+ scale_fill_manual(values=c("deepskyblue4", "firebrick2"))
p_nb_1e6

compare_means(nb_1e7 ~ tsg, data = df)
p <- ggplot(data = df, aes(y=nb_1e7,x=tsg,fill=tsg))+geom_boxplot(outlier.shape = NA)
ylim1 = boxplot.stats(df$nb_1e7)$stats[c(1, 5)]
p_nb_1e7 = p + coord_cartesian(ylim = ylim1*1.5) + stat_compare_means(label.x = 0.6, label.y = 12)+ scale_fill_manual(values=c("deepskyblue4", "firebrick2"))
p_nb_1e7


# local ranks
compare_means(local_1e6 ~ tsg, data = df)
p <- ggplot(data = df, aes(y=local_1e6,x=tsg,fill=tsg))+geom_boxplot(outlier.shape = NA)
ylim1 = boxplot.stats(df$local_1e6)$stats[c(1, 5)]
p_lc_1e6 = p + coord_cartesian(ylim = ylim1*1.5) + stat_compare_means(label.x = 1.6, label.y = 400)+ scale_fill_manual(values=c("deepskyblue4", "firebrick2"))
p_lc_1e6

compare_means(local_1e7 ~ tsg, data = df)
p <- ggplot(data = df, aes(y=local_1e7,x=tsg,fill=tsg))+geom_boxplot(outlier.shape = NA)
ylim1 = boxplot.stats(df$local_1e7)$stats[c(1, 5)]
p_lc_1e7 = p + coord_cartesian(ylim = ylim1*1.5) + stat_compare_means(label.x = 1.6, label.y = 400)+ scale_fill_manual(values=c("deepskyblue4", "firebrick2"))
p_lc_1e7

#counts

compare_means(counts_1e6 ~ tsg, data = df)
p <- ggplot(data = df, aes(y=counts_1e6,x=tsg,fill=tsg))+geom_boxplot(outlier.shape = NA)
ylim1 = boxplot.stats(df$counts_1e6)$stats[c(1, 5)]
p_counts_1e6 = p + coord_cartesian(ylim = ylim1*1.5) + stat_compare_means(label.x = 0.6, label.y = 10)+ scale_fill_manual(values=c("deepskyblue4", "firebrick2"))
p_counts_1e6

compare_means(counts_1e7 ~ tsg, data = df)
p <- ggplot(data = df, aes(y=counts_1e7,x=tsg,fill=tsg))+geom_boxplot(outlier.shape = NA)
ylim1 = boxplot.stats(df$counts_1e7)$stats[c(1, 5)]
p_counts_1e7 = p + coord_cartesian(ylim = ylim1*1.5) + stat_compare_means(label.x = 1, label.y = 100)+ scale_fill_manual(values=c("deepskyblue4", "firebrick2"))
p_counts_1e7


#peak

compare_means(peak_1e6 ~ tsg, data = df)
p <- ggplot(data = df, aes(y=peak_1e6,x=tsg,fill=tsg))+geom_boxplot(outlier.shape = NA)
ylim1 = boxplot.stats(df$peak_1e6)$stats[c(1, 5)]
p_peak_1e6 = p + coord_cartesian(ylim = ylim1*20) + stat_compare_means(label.x = 0.6, label.y = 10)+ scale_fill_manual(values=c("deepskyblue4", "firebrick2"))
p_peak_1e6

compare_means(peak_1e7 ~ tsg, data = df)
p <- ggplot(data = df, aes(y=peak_1e7,x=tsg,fill=tsg))+geom_boxplot(outlier.shape = NA)
ylim1 = boxplot.stats(df$peak_1e7)$stats[c(1, 5)]
p_peak_1e7 = p + coord_cartesian(ylim = ylim1*20) + stat_compare_means(label.x = 1, label.y = 100)+ scale_fill_manual(values=c("deepskyblue4", "firebrick2"))
p_peak_1e7


### first important output images ###
a = grid.arrange(p_fs_1e6,p_fs_1e7,p_nb_1e6,p_nb_1e7,p_lc_1e6,p_lc_1e7,p_counts_1e6,p_counts_1e7,p_peak_1e6,p_peak_1e7,ncol=4)



#dfx <- merge(df,cgc,by.x="gene",by.y="Gene Symbol",all=T)
dfx <- merge(df,cgc,by.x="gene",by.y="Gene Symbol",all.x =T)
dfx <- dfx[is.na(dfx$chr)==FALSE,]
names(dfx)[19] = "role"
names(dfx)[18] = "chr_band"
ylim1 = boxplot.stats(df$peak_1e7)$stats[c(1, 5)]
p_roles_ps_1e7 <- ggplot(data = dfx, aes(y=peak_1e7,x=role,fill=role))+geom_boxplot(outlier.shape = NA)+
  coord_cartesian(ylim = ylim1*8) + ggtitle("Peak score for genes (CNVs <= 1e7)")
  
ylim1 = boxplot.stats(df$peak_1e6)$stats[c(1, 5)]
p_roles_ps_1e6 <- ggplot(data = dfx, aes(y=peak_1e6,x=role,fill=role))+geom_boxplot(outlier.shape = NA)+
  coord_cartesian(ylim = ylim1*12) + ggtitle("Peak score for genes (CNVs <= 1e6)")

### second important output images ###
b = grid.arrange(p_roles_ps_1e6,p_roles_ps_1e7)

# gene length related
dfx$log.length <- log(dfx$end-dfx$start)
dfx$fixed.score <- dfx$fscore_1e6/dfx$log.length

compare_means(fixed.score ~ tsg, data = dfx)
k1 <- ggplot(data = dfx, aes(y=fixed.score,x=tsg,fill=tsg))+geom_boxplot(outlier.shape = NA)+
  coord_cartesian(ylim = ylim1*1.05) + ggtitle("Normalized peak score for genes (CNVs <= 1e6)") + stat_compare_means(label.x = 0.6, label.y = 1.5)

compare_means(peak_1e6 ~ tsg, data = dfx)
k2 <- ggplot(data = dfx, aes(y=peak_1e6,x=tsg,fill=tsg))+geom_boxplot(outlier.shape = NA)+
  coord_cartesian(ylim = ylim1*10) + ggtitle("Peak score for genes (CNVs <= 1e6)") + stat_compare_means(label.x = 0.6, label.y = 10)

dfx$fixed.score.1e7 <- dfx$fscore_1e7/dfx$log.length
compare_means(fixed.score.1e7 ~ tsg, data = dfx)
k71 <- ggplot(data = dfx, aes(y=fixed.score.1e7,x=tsg,fill=tsg))+geom_boxplot(outlier.shape = NA)+
  coord_cartesian(ylim = ylim1*3) + ggtitle("Normalized peak score for genes (CNVs <= 1e7)") + stat_compare_means(label.x = 0.6, label.y = 2.5)

compare_means(peak_1e7 ~ tsg, data = dfx)
k72 <- ggplot(data = dfx, aes(y=peak_1e7,x=tsg,fill=tsg))+geom_boxplot(outlier.shape = NA)+
  coord_cartesian(ylim = ylim1*15) + ggtitle("Peak score for genes (CNVs <= 1e7)") + stat_compare_means(label.x = 0.6, label.y = 20)

### third important output images ###
grid.arrange(k1,k2,k71,k72,nrow=2)



# peak score calculator

pdf("preliminary.pdf",width = 20,height = 16)

grid.arrange(p_fs_1e6,p_fs_1e7,p_nb_1e6,p_nb_1e7,p_lc_1e6,p_lc_1e7,p_counts_1e6,p_counts_1e7,p_peak_1e6,p_peak_1e7,ncol=4)
grid.arrange(p_roles_ps_1e6,p_roles_ps_1e7)
grid.arrange(k1,k2,k71,k72,nrow=2)

dev.off()


dfx$rank_fscore_1e6 <- NA
dfx$rank_fscore_1e6[order(dfx$fscore_1e6,decreasing = T)] <- 1:nrow(dfx)

dfx$rank_fscore_1e7 <- NA
dfx$rank_fscore_1e7[order(dfx$fscore_1e7,decreasing = T)] <- 1:nrow(dfx)

dfx$rank_peak_1e6 <- NA
dfx$rank_peak_1e6[order(dfx$peak_1e6,decreasing = T)] <- 1:nrow(dfx)

dfx$rank_peak_1e7 <- NA
dfx$rank_peak_1e7[order(dfx$peak_1e7,decreasing = T)] <- 1:nrow(dfx)

dfx$rank_counts_1e6 <- NA
dfx$rank_counts_1e6[order(dfx$counts_1e6,decreasing = T)] <- 1:nrow(dfx)

dfx$rank_counts_1e7 <- NA
dfx$rank_counts_1e7[order(dfx$counts_1e7,decreasing = T)] <- 1:nrow(dfx)

df.cgc <- dfx[is.na(dfx$role)==FALSE,]



library(reshape)
library(reshape2)

df.cgc.t <- df.cgc[,c("gene","role","rank_fscore_1e6","rank_fscore_1e7","rank_peak_1e6","rank_peak_1e7","rank_counts_1e6","rank_counts_1e7")]
df.cgc.t <- melt(df.cgc.t,id=c("gene","role"))
df.cgc.t <- df.cgc.t[df.cgc.t$role%like%"TSG"==TRUE,]
ggplot(data=df.cgc.t,aes(x=variable,y=value,group=gene,color=role))+geom_line()
ranks_comparison <- ggplot(data=df.cgc.t,aes(x=variable,y=value,color=variable))+
  geom_violin(aes(fill = variable)) + 
  geom_boxplot(width = 0.2)+
  ggtitle("ranks of fscore/peak_score/counts")+
  theme(plot.title = element_text(hjust = 0.5))



# peak score calculator

pdf("preliminary.pdf",width = 20,height = 16)

grid.arrange(p_fs_1e6,p_fs_1e7,p_nb_1e6,p_nb_1e7,p_lc_1e6,p_lc_1e7,p_counts_1e6,p_counts_1e7,p_peak_1e6,p_peak_1e7,ncol=4)
grid.arrange(p_roles_ps_1e6,p_roles_ps_1e7)
grid.arrange(k1,k2,k71,k72,nrow=2)
grid.arrange(ranks_comparison)

dev.off()








