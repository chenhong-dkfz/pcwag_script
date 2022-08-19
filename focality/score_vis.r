library(data.table)
library(ggplot2)

setwd("/abi/data/chen/newhome/promoter_focallity/data/")
f.loss <- fread("loss.csv",data.table = F)
#gene.list <- fread("/ibios/co02/chen/pcawg_2018/result/matrix/coding/LGG-US2.bed",data.table = F)
gene.list <- fread("/abi/data/chen/co02/pcawg_2018/result/matrix/coding/LGG-US2.bed",data.table = F)
gene.list <- gene.list[,c(1,2,3,4)]
gene.list <- cbind(gene.list,f.loss)
gene.list$all <- rowSums(f.loss)
gene.list$range <- gene.list$gene_end-gene.list$gene_start
gene.list$contribution <- gene.list$all/gene.list$range

cgc <- fread("/abi/data/chen/newhome/basic/cgc.txt",data.table = F)
gene.list$cgc <- "none"
gene.list[gene.list$gene_name%in%cgc$`Gene Symbol`,"cgc"] <- "cgc"
gene.list$tsg <- "none"
gene.list[gene.list$gene_name%in%cgc[cgc$`Role in Cancer`=="TSG","Gene Symbol"],"tsg"] <- "tsg"

gene.lite <- gene.list[,c(4,52,53,56)]



foc.plot <- ggplot(gene.lite, aes(x=all, y=range)) + geom_point(aes(col=tsg),shape=18)+
  ggtitle("focallity score vs gene length") +
  geom_smooth()

foc.box <-  ggplot(gene.lite, aes(x=tsg, y=all)) + geom_boxplot()

wilcox.test(data=gene.lite,all~tsg)

test.df <- gene.lite[gene.lite$tsg=="tsg",]

gene.list$delta <- 0
for(i in 2:nrow(gene.list)){gene.list[i,"delta"]<-(gene.list[i,"all"]-gene.list[i-1,"all"])}
gene.list$prev <- "none"
for(i in 2:nrow(gene.list)){gene.list[i,"prev"]<-(gene.list[i-1,"gene_name"])}
gene.list$prev_cgc <- "none"
gene.list[gene.list$prev%in%cgc$`Gene Symbol`,"prev_cgc"] <- "cgc"
gene.deep <- gene.list


foc.plot <- ggplot(gene.list, aes(x=delta, y=range)) + geom_point(aes(col=tsg),shape=18)+
  ggtitle("focallity score diff vs gene length") +
  geom_smooth()
foc.box <-  ggplot(gene.list, aes(x=tsg, y=delta)) + geom_boxplot()

################
# all together #
################

# all does not equal to sum of individual cohorts!!!

f.all <- fread("all_cohorts.txt",data.table = F)
gene.list <- fread("/abi/data/chen/co02/pcawg_2018/result/matrix/coding/LGG-US2.bed",data.table = F)
gene.list <- gene.list[,c(1,2,3,4)]
gene.list <- cbind(gene.list,f.all)
gene.list$range <- gene.list$gene_end-gene.list$gene_start
gene.list$contribution <- gene.list$all/gene.list$range

cgc <- fread("/abi/data/chen/newhome/basic/cgc.txt",data.table = F)
gene.list$cgc <- "none"
gene.list[gene.list$gene_name%in%cgc$`Gene Symbol`,"cgc"] <- "cgc"
gene.list$tsg <- "none"
gene.list[gene.list$gene_name%in%cgc[cgc$`Role in Cancer`=="TSG","Gene Symbol"],"tsg"] <- "tsg"
wilcox.test(data=gene.list,all~cgc)
wilcox.test(data=gene.list,all~tsg)

gene.list$delta <- 0
for(i in 2:nrow(gene.list)){gene.list[i,"delta"]<-(gene.list[i,"all"]-gene.list[i-1,"all"])}
gene.list$prev <- "none"
for(i in 2:nrow(gene.list)){gene.list[i,"prev"]<-(gene.list[i-1,"gene_name"])}
gene.list$prev_cgc <- "none"
gene.list[gene.list$prev%in%cgc$`Gene Symbol`,"prev_cgc"] <- "cgc"
gene.list$prev_chrom <- "none"
for(i in 2:nrow(gene.list)){gene.list[i,"prev_chrom"]<-(gene.list[i-1,"#chr"])}

gene.list$foll_delta <-0
for(i in 1:(nrow(gene.list)-1)){gene.list[i,"foll_delta"]<-(gene.list[i,"all"]-gene.list[i+1,"all"])}

gene.list$dev <- 0

# opt 1
gene.list$dev <- 0.5*(gene.list$delta+gene.list$foll_delta)
gene.list$edge <- "false"
for(i in 2:(nrow(gene.list)-1)){
  if(gene.list[i,"#chr"]!=gene.list[i-1,"#chr"]|gene.list[i,"#chr"]!=gene.list[i+1,"#chr"]){
    gene.list[i,"edge"] <- "true"
  }
}
gene.list[1,"edge"] <- "true"
gene.list[nrow(gene.list),"edge"] <- "true"
gene.ftd <- gene.list[gene.list$edge=="false",]

wilcox.test(data=gene.ftd,dev~cgc)
wilcox.test(data=gene.ftd,dev~tsg)
ggplot(data = gene.ftd,aes(x=cgc,y=dev))+geom_boxplot()+ylim(0,100)
vs.plot <- ggplot(data = gene.ftd,aes(x=tsg,y=dev,fill=tsg))+geom_boxplot()+ylim(0,100)+ 
            stat_compare_means(label.x = 1.5, label.y = 100)+ggtitle("average scores between TSGs and normal genes")

gene.ftd <- gene.ftd[,c(1,4,5,6,8,9,15)]
gene.ftd <- gene.ftd[gene.ftd$dev>0,]

# opt 2
for(i in 1:(nrow(gene.list))){
  if(gene.list[i,"foll_delta"]>gene.list[i,"delta"]){
    gene.list[i,"dev"]=gene.list[i,"foll_delta"]
  }else{
    gene.list[i,"dev"]=gene.list[i,"delta"]
  }
}


gene.all <- gene.list

# write.table(gene.all,"/abi/data/chen/newhome/promoter_focallity/data/results_all_genes.csv",col.names=T,row.names=F,quote=F,sep="\t")

gene.locate <- function(col.chr,col.start){
  if(col.chr=="1"){cen=125000000;fin=247200000}
  if(col.chr=="2"){cen=93300000;fin=242800000}
  if(col.chr=="3"){cen=91000000;fin=199400000}
  if(col.chr=="4"){cen=50400000;fin=191300000}
  if(col.chr=="5"){cen=48400000;fin=180800000}
  if(col.chr=="6"){cen=61000000;fin=170900000}
  if(col.chr=="7"){cen=59900000;fin=158800000}
  if(col.chr=="8"){cen=45600000;fin=146300000}
  if(col.chr=="9"){cen=49000000;fin=140400000}
  if(col.chr=="10"){cen=40200000;fin=135400000}
  if(col.chr=="11"){cen=53700000;fin=134500000}
  if(col.chr=="12"){cen=35800000;fin=132300000}
  if(col.chr=="13"){cen=17900000;fin=114100000}
  if(col.chr=="14"){cen=17600000;fin=106300000}
  if(col.chr=="15"){cen=19000000;fin=100300000}
  if(col.chr=="16"){cen=36600000;fin=88800000}
  if(col.chr=="17"){cen=24000000;fin=78700000}
  if(col.chr=="18"){cen=17200000;fin=76100000}
  if(col.chr=="19"){cen=26500000;fin=63800000}
  if(col.chr=="20"){cen=27500000;fin=62400000}
  if(col.chr=="21"){cen=13200000;fin=46900000}
  if(col.chr=="22"){cen=14700000;fin=49500000}
  if(col.chr=="X"){cen=60600000;fin=154900000}
  if(col.chr=="Y"){cen=12500000;fin=57700000}
  if(col.start<cen){
    loc <- col.start/cen
    arms <- "p"
  }else{
    loc <- (col.start-cen)/(fin-cen)
    arms <- "q"
  }
  return(c(loc,arms))
}

gene.all$loc <- 0
gene.all$arm <- "na"

for(i in 1:nrow(gene.all)){
chrin <- gene.all[i,1]
sttin <- gene.all[i,2]
result <- gene.locate(chrin,sttin)
gene.all[i,"loc"] <- result[1]
gene.all[i,"arm"] <- result[2]
}

# rearrange to gene.all
colnames(gene.all)[1] <- "chr"
colnames(gene.all)[5] <- "score"
colnames(gene.all)[6] <- "length"

gene.all.new <- gene.all[,c(1,2,3,4,5,6,8,9,16,17,18)]
x0 <- gene.all.new[,1:6]
y0 <- x0[c(1,1:(nrow(x0)-1)),]
z0 <- x0[c(2:nrow(x0),nrow(x0)),]
df.head <- cbind(x0,y0,z0)
df.head <- df.head[,c(1,2,3,4,5,6,7,10,11,12,13,16,17,18)]
colnames(df.head) <- c("chr","start","end","gene","score","length",
                      "prev_chr","prev_gene","prev_score","prev_length",
                      "next_chr","next_gene","next_score","next_length")
df.tail <- gene.all.new[,c(7,8,9,10,11)]
df <- cbind(df.head,df.tail)
#filter edges
df <- df[df$edge!="true",]
df <- df[,c(1,2,3,4,5,6,8,9,10,12,13,14,15,16,18,19)]
df$delta_1 <- df$score-df$prev_score
df$delta_2 <- df$score-df$next_score
df$peak <- "false"
df[df$delta_1>0&df$delta_2>0,"peak"] <- "true"
df$peak_value <- pmax(df$delta_1,df$delta_2)

df$rel_score <- df$score/log10(df$length)
df$rel_prev_score <- df$prev_score/log10(df$prev_length)
df$rel_next_score <- df$next_score/log10(df$next_length)
df$rel_peak <- "false"
df$rel_delta_1 <- df$rel_score-df$rel_prev_score
df$rel_delta_2 <- df$rel_score-df$rel_next_score

df[df$rel_delta_1>0&df$rel_delta_2>0,"rel_peak"] <- "true"
df$rel_peak_value <-pmax(df$rel_delta_1,df$rel_delta_2)

write.table(df,"/abi/data/chen/newhome/promoter_focallity/result/results_31042019.csv",col.names=T,row.names=F,quote=F,sep="\t")
df.lite <- df[,c(1,4,5,6,7,10,13,14,15,19,20,21,24,27)]
write.table(df.lite,"/abi/data/chen/newhome/promoter_focallity/result/results_lite_31042019.csv",col.names=T,row.names=F,quote=F,sep="\t")

#write.table(gene.all,"/abi/data/chen/newhome/promoter_focallity/result/results_all_genes.csv",col.names=T,row.names=F,quote=F,sep="\t")

library(ggrepel,lib="/abi/data/chen/R/lib64/3.5.1")

df.stat <- df.lite[,c(1,2,3,4,7,8,9,11,12,14)]
df.stat$loc <- as.numeric(as.character(df.stat$loc))
sc <- ggplot(data=df.stat,aes(y=score,x=loc))+geom_point(fill=chr)
rel_vs_loc <- ggplot(data=df.stat,aes(y=rel_score,x=loc,color=chr))+geom_point()
sc <- ggplot(data=df.stat,aes(y=rel_peak_value,x=loc,color=chr,
                              label=ifelse( (rel_peak_value>50&loc>0.01&loc<0.99)  ,as.character(gene),'')))+geom_point() + geom_label_repel() 
  
geom_text(aes(label=ifelse( (rel_peak_value>50&loc>0.01&loc<0.99)  ,as.character(gene),'')),hjust=0,vjust=0)

##############
# quick load #
##############
g.lite <- fread("/abi/data/chen/newhome/promoter_focallity/result/results_lite_31042019.csv",data.table=F)
plot0 <- ggplot(data=g.lite,aes(y=rel_peak_value,x=loc,color=chr,
                                label=ifelse( (rel_peak_value>50&loc>0.01&loc<0.99)  ,as.character(gene),'')))+
                                geom_point() + geom_label_repel() 
  
  
  
  


com.list <- fread("raw.txt",data.table = F)
gene.com <- fread("/ibios/co02/chen/pcawg_2018/result/matrix/coding/LGG-US2.bed",data.table = F)
gene.com <- gene.com[,c(1,2,3,4)]
gene.com <- cbind(gene.com,com.list)



###############
# new version #
###############

# all does not equal to sum of individual cohorts!!!

f.all <- fread("all_cohorts_new.txt",data.table = F)
gene.list <- fread("/abi/data/chen/co02/pcawg_2018/result/matrix/coding/LGG-US2.bed",data.table = F)
gene.list <- gene.list[,c(1,2,3,4)]
gene.list <- cbind(gene.list,f.all)
gene.list$range <- gene.list$gene_end-gene.list$gene_start
gene.list$contribution <- gene.list$all/gene.list$range

cgc <- fread("/abi/data/chen/newhome/basic/cgc.txt",data.table = F)
gene.list$cgc <- "none"
gene.list[gene.list$gene_name%in%cgc$`Gene Symbol`,"cgc"] <- "cgc"
gene.list$tsg <- "none"
gene.list[gene.list$gene_name%in%cgc[cgc$`Role in Cancer`=="TSG","Gene Symbol"],"tsg"] <- "tsg"
wilcox.test(data=gene.list,all~cgc)
wilcox.test(data=gene.list,all~tsg)
names(gene.list)[1] <- "chr"

gene.list$edge <- "none"
for(i in 2:nrow(gene.list)){
  if(gene.list[i,"chr"]!=gene.list[i-1,"chr"]){
    gene.list[i,"edge"] = "head"
  }
}
for(i in 1:(nrow(gene.list)-1)){
  if(gene.list[i,"chr"]!=gene.list[i+1,"chr"]){
    gene.list[i,"edge"] = "tail"
  }
}
gene.list[1,"edge"] <- "head"
gene.list[nrow(gene.list),"edge"] <- "tail"

gene.list$edge.score <- 0
for (i in 1:nrow(gene.list)){
  if(gene.list[i,"edge"]=="none"){
    gene.list[i,"edge.score"]=0.5*(2*gene.list[i,"all"]-gene.list[i-1,"all"]-gene.list[i+1,"all"])
  }else if(gene.list[i,"edge"]=="head"){
    gene.list[i,"edge.score"]=0.5*(2*gene.list[i,"all"]-gene.list[i+1,"all"]-gene.list[i+1,"all"])
  }else if(gene.list[i,"edge"]=="tail"){
    gene.list[i,"edge.score"]=0.5*(2*gene.list[i,"all"]-gene.list[i-1,"all"]-gene.list[i-1,"all"])
  }
}


gene.list$rk.range <- "exclusive"
for (i in 6:(nrow(gene.list)-5)){
  if(gene.list[i,"chr"]==gene.list[i+5,"chr"] & gene.list[i,"chr"]==gene.list[i-5,"chr"]){
    gene.list[i,"rk.range"]="inclusive"
  }
}

gene.list$rank <- 0
for (i in 1:nrow(gene.list)){
  if(gene.list[i,"rk.range"]=="inclusive"){
    temp_list = c(gene.list[i-5,"all"],gene.list[i-4,"all"],gene.list[i-3,"all"],gene.list[i-2,"all"],gene.list[i-1,"all"],
                  gene.list[i+1,"all"],gene.list[i+2,"all"],gene.list[i+3,"all"],gene.list[i+4,"all"],gene.list[i+5,"all"],
                  gene.list[i,"all"])
    temp_list = sort(temp_list,decreasing = TRUE)
    gene.list[i,"rank"] = match(gene.list[i,"all"],temp_list)
  }
}

gene.list$fs <- "none"
fs.list <- fread("/abi/data/chen/newhome/fragile_site_bed/gene.txt",data.table = F)
names(fs.list) <- c("chr","start","end","gene","score","strand")
gene.list[gene.list$gene_name%in%fs.list$`gene`,"fs"] <- "fragile"

# visualization
ggplot(data = gene.list,aes(x=all,y=edge.score,alpha=tsg,shape=fs))+geom_point()
ggplot(gene.list, aes(x=tsg, y=all,color=tsg)) + geom_boxplot()

ylim1 = boxplot.stats(gene.list$edge.score)$stats[c(1, 5)]
p0 <- ggplot(gene.list, aes(x=cgc, y=edge.score,color=cgc)) + geom_boxplot()
p0 + coord_cartesian(ylim = ylim1*20)

ggplot(gene.list, aes(x=cgc, y=rank,color=cgc)) + geom_boxplot()

write.table(gene.list,file="focallity_scores_threshold_10e6.tsv",sep="\t",quote = F,row.names = F)


gene.list$delta <- 0
for(i in 2:nrow(gene.list)){gene.list[i,"delta"]<-(gene.list[i,"all"]-gene.list[i-1,"all"])}
gene.list$prev <- "none"
for(i in 2:nrow(gene.list)){gene.list[i,"prev"]<-(gene.list[i-1,"gene_name"])}
gene.list$prev_cgc <- "none"
gene.list[gene.list$prev%in%cgc$`Gene Symbol`,"prev_cgc"] <- "cgc"
gene.list$prev_chrom <- "none"
for(i in 2:nrow(gene.list)){gene.list[i,"prev_chrom"]<-(gene.list[i-1,"#chr"])}

gene.list$foll_delta <-0
for(i in 1:(nrow(gene.list)-1)){gene.list[i,"foll_delta"]<-(gene.list[i,"all"]-gene.list[i+1,"all"])}

gene.list$dev <- 0

# opt 1
gene.list$dev <- 0.5*(gene.list$delta+gene.list$foll_delta)
gene.list$edge <- "false"
for(i in 2:(nrow(gene.list)-1)){
  if(gene.list[i,"#chr"]!=gene.list[i-1,"#chr"]|gene.list[i,"#chr"]!=gene.list[i+1,"#chr"]){
    gene.list[i,"edge"] <- "true"
  }
}
gene.list[1,"edge"] <- "true"
gene.list[nrow(gene.list),"edge"] <- "true"
gene.ftd <- gene.list[gene.list$edge=="false",]

wilcox.test(data=gene.ftd,dev~cgc)
wilcox.test(data=gene.ftd,dev~tsg)
ggplot(data = gene.ftd,aes(x=cgc,y=dev))+geom_boxplot()+ylim(0,100)
vs.plot <- ggplot(data = gene.ftd,aes(x=tsg,y=dev,fill=tsg))+geom_boxplot()+ylim(0,100)+ 
  stat_compare_means(label.x = 1.5, label.y = 100)+ggtitle("average scores between TSGs and normal genes")

gene.ftd <- gene.ftd[,c(1,4,5,6,8,9,15)]
gene.ftd <- gene.ftd[gene.ftd$dev>0,]

# opt 2
for(i in 1:(nrow(gene.list))){
  if(gene.list[i,"foll_delta"]>gene.list[i,"delta"]){
    gene.list[i,"dev"]=gene.list[i,"foll_delta"]
  }else{
    gene.list[i,"dev"]=gene.list[i,"delta"]
  }
}


gene.all <- gene.list

# write.table(gene.all,"/abi/data/chen/newhome/promoter_focallity/data/results_all_genes.csv",col.names=T,row.names=F,quote=F,sep="\t")

gene.locate <- function(col.chr,col.start){
  if(col.chr=="1"){cen=125000000;fin=247200000}
  if(col.chr=="2"){cen=93300000;fin=242800000}
  if(col.chr=="3"){cen=91000000;fin=199400000}
  if(col.chr=="4"){cen=50400000;fin=191300000}
  if(col.chr=="5"){cen=48400000;fin=180800000}
  if(col.chr=="6"){cen=61000000;fin=170900000}
  if(col.chr=="7"){cen=59900000;fin=158800000}
  if(col.chr=="8"){cen=45600000;fin=146300000}
  if(col.chr=="9"){cen=49000000;fin=140400000}
  if(col.chr=="10"){cen=40200000;fin=135400000}
  if(col.chr=="11"){cen=53700000;fin=134500000}
  if(col.chr=="12"){cen=35800000;fin=132300000}
  if(col.chr=="13"){cen=17900000;fin=114100000}
  if(col.chr=="14"){cen=17600000;fin=106300000}
  if(col.chr=="15"){cen=19000000;fin=100300000}
  if(col.chr=="16"){cen=36600000;fin=88800000}
  if(col.chr=="17"){cen=24000000;fin=78700000}
  if(col.chr=="18"){cen=17200000;fin=76100000}
  if(col.chr=="19"){cen=26500000;fin=63800000}
  if(col.chr=="20"){cen=27500000;fin=62400000}
  if(col.chr=="21"){cen=13200000;fin=46900000}
  if(col.chr=="22"){cen=14700000;fin=49500000}
  if(col.chr=="X"){cen=60600000;fin=154900000}
  if(col.chr=="Y"){cen=12500000;fin=57700000}
  if(col.start<cen){
    loc <- col.start/cen
    arms <- "p"
  }else{
    loc <- (col.start-cen)/(fin-cen)
    arms <- "q"
  }
  return(c(loc,arms))
}

gene.all$loc <- 0
gene.all$arm <- "na"

for(i in 1:nrow(gene.all)){
  chrin <- gene.all[i,1]
  sttin <- gene.all[i,2]
  result <- gene.locate(chrin,sttin)
  gene.all[i,"loc"] <- result[1]
  gene.all[i,"arm"] <- result[2]
}

# rearrange to gene.all
colnames(gene.all)[1] <- "chr"
colnames(gene.all)[5] <- "score"
colnames(gene.all)[6] <- "length"

gene.all.new <- gene.all[,c(1,2,3,4,5,6,8,9,16,17,18)]
x0 <- gene.all.new[,1:6]
y0 <- x0[c(1,1:(nrow(x0)-1)),]
z0 <- x0[c(2:nrow(x0),nrow(x0)),]
df.head <- cbind(x0,y0,z0)
df.head <- df.head[,c(1,2,3,4,5,6,7,10,11,12,13,16,17,18)]
colnames(df.head) <- c("chr","start","end","gene","score","length",
                       "prev_chr","prev_gene","prev_score","prev_length",
                       "next_chr","next_gene","next_score","next_length")
df.tail <- gene.all.new[,c(7,8,9,10,11)]
df <- cbind(df.head,df.tail)
#filter edges
df <- df[df$edge!="true",]
df <- df[,c(1,2,3,4,5,6,8,9,10,12,13,14,15,16,18,19)]
df$delta_1 <- df$score-df$prev_score
df$delta_2 <- df$score-df$next_score
df$peak <- "false"
df[df$delta_1>0&df$delta_2>0,"peak"] <- "true"
df$peak_value <- pmax(df$delta_1,df$delta_2)

df$rel_score <- df$score/log10(df$length)
df$rel_prev_score <- df$prev_score/log10(df$prev_length)
df$rel_next_score <- df$next_score/log10(df$next_length)
df$rel_peak <- "false"
df$rel_delta_1 <- df$rel_score-df$rel_prev_score
df$rel_delta_2 <- df$rel_score-df$rel_next_score

df[df$rel_delta_1>0&df$rel_delta_2>0,"rel_peak"] <- "true"
df$rel_peak_value <-pmax(df$rel_delta_1,df$rel_delta_2)

write.table(df,"/abi/data/chen/newhome/promoter_focallity/result/results_31042019.csv",col.names=T,row.names=F,quote=F,sep="\t")
df.lite <- df[,c(1,4,5,6,7,10,13,14,15,19,20,21,24,27)]
write.table(df.lite,"/abi/data/chen/newhome/promoter_focallity/result/results_lite_31042019.csv",col.names=T,row.names=F,quote=F,sep="\t")

#write.table(gene.all,"/abi/data/chen/newhome/promoter_focallity/result/results_all_genes.csv",col.names=T,row.names=F,quote=F,sep="\t")

library(ggrepel,lib="/abi/data/chen/R/lib64/3.5.1")

df.stat <- df.lite[,c(1,2,3,4,7,8,9,11,12,14)]
df.stat$loc <- as.numeric(as.character(df.stat$loc))
sc <- ggplot(data=df.stat,aes(y=score,x=loc))+geom_point(fill=chr)
rel_vs_loc <- ggplot(data=df.stat,aes(y=rel_score,x=loc,color=chr))+geom_point()
sc <- ggplot(data=df.stat,aes(y=rel_peak_value,x=loc,color=chr,
                              label=ifelse( (rel_peak_value>50&loc>0.01&loc<0.99)  ,as.character(gene),'')))+geom_point() + geom_label_repel() 

geom_text(aes(label=ifelse( (rel_peak_value>50&loc>0.01&loc<0.99)  ,as.character(gene),'')),hjust=0,vjust=0)

##############
# quick load #
##############
g.lite <- fread("/abi/data/chen/newhome/promoter_focallity/result/results_lite_31042019.csv",data.table=F)
plot0 <- ggplot(data=g.lite,aes(y=rel_peak_value,x=loc,color=chr,
                                label=ifelse( (rel_peak_value>50&loc>0.01&loc<0.99)  ,as.character(gene),'')))+
  geom_point() + geom_label_repel() 






com.list <- fread("raw.txt",data.table = F)
gene.com <- fread("/ibios/co02/chen/pcawg_2018/result/matrix/coding/LGG-US2.bed",data.table = F)
gene.com <- gene.com[,c(1,2,3,4)]
gene.com <- cbind(gene.com,com.list)
