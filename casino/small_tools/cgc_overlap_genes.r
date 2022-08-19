cgc<-read.table(paste0("/ibios/co02/chen/casino/gencode/results/analysis/",
                       "cgc_list.txt"),sep="\t",header=T)
colnames(cgc)<-"gene"
casino<-read.table(paste0("/ibios/co02/chen/casino/gencode/results/analysis/",
                        "score_ranked2.bed"),sep="\t",header=F)
colnames(casino)<-c("chromosome","start","end","gene","score") 
casino$ID<-seq.int(nrow(casino))
overlap <- merge(x=cgc,y=casino,by.x="gene",by.y="gene")

mutsig <- read.table(paste0("/ibios/co02/chen/casino/gencode/results/analysis/",
                        "mutsig.tsv"),sep="\t",header=T,comment="")
colnames(mutsig)[1]<-"gene"
music <- read.table(paste0("/ibios/co02/chen/casino/gencode/results/analysis/",
                            "significantly_mutated_genes.smg"),sep="\t",header=T,comment="")
colnames(music)[1]<-"gene"

overlap_mutsig<-merge(x=mutsig,y=casino,by="gene")
overlap_music<-merge(x=music,y=casino,by="gene")

library("ggplot2")
casino2<-casino[,c(4,5,6)]
casino2$mutsig<-(rep(0,dim(casino2)[1]))
casino2$music<-(rep(0,dim(casino2)[1]))
casino2$cgc<-(rep(0,dim(casino2)[1]))
casino2[which(casino2$gene %in% overlap$gene),]$cgc = 1
casino2[which(casino2$gene %in% overlap_mutsig$gene),]$mutsig = 2
casino2[which(casino2$gene %in% overlap_music$gene),]$music = 4
casino2$status<-casino2$mutsig + casino2$music + casino2$cgc
#ggplot2.barplot(data=casino2, xName="score",groupName="status")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
df<-casino2

pdf(paste0("/ibios/co02/chen/casino/gencode/results/analysis/",
           "result_overlap_top_1000.pdf"))

ggplot(df,aes(x=ID,y=score,fill=factor(status)))+
  geom_bar(stat='identity')+
  scale_fill_manual(
    name = "result overlaps",
    breaks=c(0,1,2,3,4,5,6,7),
    labels=c("no match","cgc","mutsig","cgc&mutsig","music","cgc&music","mutsig&music","cgc&mutsig&music"),
    values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
    scale_x_discrete(breaks=df$ID,labels=df$gene)+
    theme(axis.text.x=element_text(angle=90,hjust=0,size=1))+
    labs(x="genes",y="score")

dev.off()