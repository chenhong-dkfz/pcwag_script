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

# read files

#snv.table <- fread("/home/hongc/PhD/TumorPrint/single_gene/test/snv.txt",data.table = F)
#cnv.table <- fread("/home/hongc/PhD/TumorPrint/single_gene/test/cnv.txt",data.table = F)
#exp.table <- fread("/home/hongc/PhD/TumorPrint/single_gene/test/exp.txt",data.table = F)


# input file I and II

snv.table <- fread("/ibios/co02/chen/pcawg_2018/result/file_cache/mut.csv",data.table = F)
exp.table <- fread("/ibios/co02/chen/pcawg_2018/result/file_cache/exp_file_ta.csv",data.table = F)

rownames(snv.table) <- snv.table[, 1] 
snv.table <- snv.table[, -1]  
#rownames(cnv.table) <- cnv.table[, 1] 
#cnv.table <- cnv.table[, -1]  
rownames(exp.table) <- exp.table[, 1] 
exp.table <- exp.table[, -1]  

####################
# definition tables#
####################

#input file III
pid.cohort <- fread("/ibios/co02/chen/pcawg_2018/pid_list/pid.tsv",data.table = F,header = F)
rownames(pid.cohort) <- pid.cohort[, 1] 
#pid.cohort <- data.frame(c("P1","P2","P3","P4"),c("C1","C1","C2","C2"))
colnames(pid.cohort) <- c("pid","cohort")
#rownames(pid.cohort) <- pid.cohort[,1]

# gene names

args = commandArgs(trailingOnly=TRUE)
gene.name <- args[1]
output.dir <- args[2]
#gene.name <- "REV1"

# combine tables

gene.matrix  <- NULL

gene.matrix <- rbind(gene.matrix ,snv.table[gene.name,])
#gene.matrix <- rbind(gene.matrix ,cnv.table[gene.name,])
gene.matrix <- smartbind(gene.matrix ,exp.table[gene.name,])
#rownames(gene.matrix) <- c("snv","cnv","exp")
rownames(gene.matrix) <- c("snv","exp")

mutation.split<- function(mutation,annotation){
  gene.matrix <<- rbind(gene.matrix,rep(0,ncol(gene.matrix)))
  rownames(gene.matrix)[nrow(gene.matrix)] <<- annotation 
  add_line <- grepl(annotation,gene.matrix[mutation,])
  add_line <- as.numeric(add_line) 
  gene.matrix[annotation,] <<- add_line
}

mutation.biallelic<- function(){
  gene.matrix <<- rbind(gene.matrix,rep(0,ncol(gene.matrix)))
  rownames(gene.matrix)[nrow(gene.matrix)] <<- "biallelic"
  gene.matrix <<- rbind(gene.matrix,rep(0,ncol(gene.matrix)))
  rownames(gene.matrix)[nrow(gene.matrix)] <<- "monoallelic"
  gene.matrix <<- rbind(gene.matrix,rep(0,ncol(gene.matrix)))
  rownames(gene.matrix)[nrow(gene.matrix)] <<- "cnv gain"
  
  add_line_snv <- grepl("SNV_4|SNV_5|INDEL_2|INDEL_3",gene.matrix["snv",]) # have snv
  add_line_cnd <- grepl("CNV_LOSS_HET",gene.matrix["snv",]) # have cnv loss
  add_line_cng <- grepl("CNV_GAIN",gene.matrix["snv",]) # have cnv gain
  add_line_hom <- grepl("CNV_LOSS_HOM",gene.matrix["snv",]) # have cnv loss hom
  
  add_line_bad <- (add_line_snv & add_line_cnd) | add_line_hom
  add_line_mad <- (add_line_snv | add_line_cnd) & (!add_line_bad)
  
  add_line_bad <- as.numeric(add_line_bad)
  add_line_mad <- as.numeric(add_line_mad)
  add_line_cng <- as.numeric(add_line_cng)

  gene.matrix["biallelic",] <<- add_line_bad
  gene.matrix["monoallelic",] <<- add_line_mad
  gene.matrix["cnv gain",] <<- add_line_cng
  
}


mutation.split("snv","CNV_GAIN")
mutation.split("snv","CNV_LOSS_HET")
mutation.split("snv","CNV_LOSS_HOM")
#mutation.split("snv","SNV_1")
#mutation.split("snv","SNV_2")
#mutation.split("snv","SNV_3")
mutation.split("snv","SNV_4")
mutation.split("snv","SNV_5")
mutation.split("snv","INDEL_2")
mutation.split("snv","INDEL_3")

mutation.biallelic()



# get cohort information
vec <- colnames(gene.matrix)
coh.line <- as.character(pid.cohort[vec,2])
gene.matrix <- rbind(gene.matrix,coh.line)
rownames(gene.matrix)[nrow(gene.matrix)] <- "cohort"

gene.matrix <- gene.matrix[!rownames(gene.matrix) %in% c("snv"),]
rownames(gene.matrix)[5] <- "nonsynonymous snv"
rownames(gene.matrix)[6] <- "Stopgain SNV"
rownames(gene.matrix)[7] <- "nonframeshift INDEL"
rownames(gene.matrix)[8] <- "frameshift INDEL"


#colour preparation
makeColors <- function(){
  maxColors <- 50
  usedColors <- c()
  possibleColors <- colorRampPalette( brewer.pal( 9 , "Set1" ) )(maxColors)
  
  function(values){
    newKeys <- setdiff(values, names(usedColors))
    newColors <- possibleColors[1:length(newKeys)]
    usedColors.new <-  c(usedColors, newColors)
    names(usedColors.new) <- c(names(usedColors), newKeys)
    usedColors <<- usedColors.new
    
    possibleColors <<- possibleColors[length(newKeys)+1:maxColors]
    usedColors
  }
} 

mkColor <- makeColors()

###################################
# remove patient without mutation #
###################################

filter.matrix <- function(mtx,filter_type){
  if(filter_type == "gain"){
    mtx <- mtx[which(mtx$CNV_GAIN!=0),]
  }
  else if(filter_type == "biallelic"){
    mtx <- mtx[which(mtx$biallelic!=0),]
  }
  else if(filter_type == "monoallelic"){
    mtx <- mtx[which(mtx$monoallelic!=0),]
  }
  else if(filter_type == "loss"){
    mtx <- mtx[which(mtx$monoallelic!=0|mtx$biallelic!=0),]
  }
  return(mtx)
}



#################
# plot funcitons#
#################

expprint <- function(exp.matrix){
  # a exp matrix has rownames of pids, two columnsL exp and cohort
  expplot <- ggplot(exp.matrix,aes(x=factor(rownames(exp.matrix),levels=rownames(exp.matrix)),y=as.numeric(as.character(exp)))) + 
    geom_bar(stat="identity",aes(fill="red")) +
    geom_hline(yintercept=0,colour="red")+
    ylab("Expression(Z-Scores)") + 
    guides(fill=guide_legend(nrow=2)) +
    theme(legend.key.width = unit(2,"mm"),
          axis.text = element_text(size=10),
          axis.title.x =  element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_line(colour = "ivory"),
          panel.background =  element_rect(fill="transparent",color = NA),
          legend.position = "none",
          legend.title = element_blank(),
          plot.margin= unit(c(0, 1, 1, 1), "lines"))
 return(expplot)
}

mutprint <- function(mut.matrix){
  # mut.part with 5 columsn:"patient","variable":mutations,"value":0 or 1
  # "cohort" and "mutation":same as colnames
  mutplot <- ggplot(mut.matrix, aes(patient,variable)) + xlab("Samples") + ylab("Alterations")+
    geom_tile(aes(col="white",fill = factor((cohort)),alpha=factor(value))) +
    scale_fill_manual(values = mkColor(factor(mut.part$cohort)))+
    scale_alpha_discrete(guide=FALSE)+
    scale_color_discrete(guide=FALSE)+
    guides(fill=guide_legend(nrow=2)) +
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_line(colour = "ivory"),
          panel.background =  element_rect(fill="white"),
          legend.position = "bottom",
          legend.title = element_blank(),
          plot.margin= unit(c(-1, 1, 0, 1), "lines"))
  return(mutplot)
}

mutprint.meta <- function(mut.matrix){
  # mut.part with 5 columsn:"patient","variable":mutations,"value":0 or 1
  # "cohort" and "mutation":same as colnames
  mutplot <- ggplot(mut.matrix, aes(patient,variable)) + xlab(cohort) + ylab("Alterations")+
    geom_tile(aes(col="white",fill = factor((project)),alpha=factor(value))) +
    #scale_fill_manual(values = mkColor(factor(cohort)))+
    scale_alpha_discrete(guide=FALSE)+
    scale_color_discrete(guide=FALSE)+
    guides(fill=guide_legend(nrow=2)) +
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_line(colour = "ivory"),
          panel.background =  element_rect(fill="white"),
          legend.position = "bottom",
          legend.title = element_blank(),
          plot.margin= unit(c(-1, 1, 0, 1), "lines"))
  return(mutplot)
}

blankprint <- function(title){
  bplot <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank() + 
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = 14,hjust = 0.5))+
    ggtitle(paste(gene.name," ",title,sep=""))
  
  return(bplot)
}

pieplot <- function(distribution,event,entropy){
  # distribution df: 2 columns:cohort and event
  pplot <- ggplot(data=distribution, aes_string(x=factor(1),y=event,fill = factor(distribution$cohort))) + 
    geom_bar(stat = "identity")+ coord_polar(theta="y")+
    ggtitle(paste(gene.name,event,"entropy=",entropy,sep="::"))+
    scale_fill_manual(values = mkColor(factor(distribution$cohort))) +   
    theme(axis.title = element_blank(),
          plot.title = element_text(size = 14,hjust = 0.5),
          panel.grid.major = element_line(colour = "grey"), 
          panel.grid.minor = element_line(colour = "grey"),
          panel.background = element_blank(),
          axis.text.x = element_text(size=20,angle = 0),
          axis.text.y = element_blank(),
          legend.title = element_blank(),
          #legend.position = "bottom",
          axis.ticks = element_blank())
  return(pplot)
}


################
# stat function#
################

generate.dict <- function(data){ # input should 
  t1 <- table(data[,2])
  t2 <- data.frame(t1)
  names(t2) <- c("cohort","size")
  return(t2)
  }



####################
# filter white list#
####################

#input file IV
list.table <- fread("/ibios/co02/chen/pcawg_2018/regions/pcawg_sample_info_v1.txt",
                    data.table = F)
coh.info <- list.table[,c(1,3)]
list.table <- list.table[,c(1,7)]
list.table <- list.table[list.table$donor_wgs_included_excluded=="Included",]
gene.matrix <- gene.matrix[,colnames(gene.matrix)%in% list.table$tumor_aliquot_id_wgs]

################
##### GAIN #####
################

df0 <- data.frame(t(gene.matrix),check.names = F)
df0 <- filter.matrix(df0,"gain")                # show type of mutation
df0 <- df0[order(df0$cohort,as.numeric(as.character(df0$exp))),] 
df0 <- as.data.frame(t(df0),check.names = F)


# make plot for expression
if (ncol(df0)!=0){
exp.part <- data.frame(t(df0),check.names = F)
exp.part <- exp.part[,c("exp","cohort")] 
exp.part$exp <- as.numeric(as.character(exp.part$exp))

q <- expprint(exp.part)

# make plot for mutations
mut.part <- data.frame(t(df0),check.names = F)
mut.part <- mut.part[,!colnames(mut.part)%in% c("exp","cohort","cnv gain")]
mut.part$patient <- rownames(mut.part)
mut.part <- melt(mut.part,id.vars="patient")
colnames(pid.cohort) <- c("patient","cohort")
mut.part <- join(mut.part,pid.cohort)
mut.part$mutation <- rownames(mut.part)
mut.part$patient <- factor(mut.part$patient, unique(as.character(mut.part$patient)))

p <- mutprint(mut.part)

# generating TUMORPrint
p1 <- ggplot_gtable(ggplot_build(p))
p2 <- ggplot_gtable(ggplot_build(q))
maxWidth = unit.pmax(p1$widths[2:5], p2$widths[2:5])
p1$widths[2:5] <- maxWidth
p2$widths[2:5] <- maxWidth
plot0 <- arrangeGrob(p2,p1, heights = c(2,3),
             top=textGrob(paste(gene.name," CNV gain and expression",sep=""), gp=gpar(fontsize=14,font=8)))
}else{
  plot0 <- blankprint("CNV gain and expression")
}


###################################
##### BIALLELIC INACTIVATIONS #####
###################################

df0 <- data.frame(t(gene.matrix),check.names = F)
df0 <- filter.matrix(df0,"biallelic")                # show type of mutation
df0 <- df0[order(df0$cohort,as.numeric(as.character(df0$exp))),] 
df0 <- as.data.frame(t(df0),check.names = F)
rownames(df0) <- c("exp","cnv gain > 8","cnv loss HET","cnv loss HOM","nonsynonymous snv",
                   "stopgain snv","nonframeshift indel","frameshift indel","biallelic inactivation",
                   "monoallelic inactivation","cnv gain","cohort")

if (ncol(df0)!=0){
# make plot for expression

exp.part <- data.frame(t(df0),check.names = F)
exp.part <- exp.part[,c("exp","cohort")] 
exp.part$exp <- as.numeric(as.character(exp.part$exp))

q <-  expprint(exp.part)

# make plot for mutations


mut.part <- data.frame(t(df0),check.names = F)
mut.part <- mut.part[,!colnames(mut.part)%in% c("exp","cohort","cnv gain")]
mut.part$patient <- rownames(mut.part)
mut.part <- melt(mut.part,id.vars="patient")
colnames(pid.cohort) <- c("patient","cohort")
mut.part <- join(mut.part,pid.cohort)
mut.part$mutation <- rownames(mut.part)
mut.part$patient <- factor(mut.part$patient, unique(as.character(mut.part$patient)))

p <- mutprint(mut.part)

# generating TUMORPrint
p1 <- ggplot_gtable(ggplot_build(p))
p2 <- ggplot_gtable(ggplot_build(q))
maxWidth = unit.pmax(p1$widths[2:5], p2$widths[2:5])
p1$widths[2:5] <- maxWidth
p2$widths[2:5] <- maxWidth
plot1 <- arrangeGrob(p2,p1, heights = c(2,3),
                      top=textGrob(paste(gene.name," inactivation and expression",sep=""), gp=gpar(fontsize=14,font=8)))
}else{
  plot1 <- blankprint("inactivation and expression")
}

################
# INACTIVATIONS#
################

df0 <- data.frame(t(gene.matrix),check.names = F)
#df0 <- filter.matrix(df0,"biallelic")
#df0 <- filter.matrix(df0,"monoallelic")
df0 <- filter.matrix(df0,"loss")# show type of mutation
dim(df0)
df0 <- df0[order(df0$cohort,as.numeric(as.character(df0$exp))),] 
df0 <- as.data.frame(t(df0),check.names = F)

inact_flag <- 0




# make plot for expression
if (ncol(df0)!=0){
  exp.part <- data.frame(t(df0),check.names = F)
  exp.part <- exp.part[,c("exp","cohort")] 
  exp.part$exp <- as.numeric(as.character(exp.part$exp))
  
  q <- expprint(exp.part)
  
  # make plot for mutations
  mut.part <- data.frame(t(df0),check.names = F)
  mut.part <- mut.part[,!colnames(mut.part)%in% c("exp","cohort","cnv gain")]
  mut.part$patient <- rownames(mut.part)
  mut.part <- melt(mut.part,id.vars="patient")
  colnames(pid.cohort) <- c("patient","cohort")
  mut.part <- join(mut.part,pid.cohort)
  mut.part$mutation <- rownames(mut.part)
  mut.part$patient <- factor(mut.part$patient, unique(as.character(mut.part$patient)))
  
  p <- mutprint(mut.part)
  
  # for this step, a general overview is already generated 
  # but no output
  
  cohorts <- unique(exp.part$cohort)
  
  # input file V
  metas <- fread("/ibios/co02/chen/pcawg_2018/meta.txt",data.table = F,header = F)
  colnames(metas) <- c("coh","met")
  
  exp.part.d <- merge(exp.part,metas,by.x = "cohort",by.y = "coh")
  cohorts <- unique(exp.part.d$met)
  colnames(exp.part.d) <- c("project","exp","cohort") 
  rownames(exp.part.d) <- rownames(exp.part)
  
  mut.part.d <- merge(mut.part,metas,by.x = "cohort",by.y = "coh")
  colnames(mut.part.d)[1] <- "project"
  colnames(mut.part.d)[6] <- "cohort"
  
  myplots <- list()
  i = 1
  for (cohort in cohorts){
    local({
      exp.part.ch <- exp.part.d[which(exp.part.d$cohort==cohort),]
      mut.part.ch <- mut.part.d[which(mut.part.d$cohort==cohort),]
      mut.part.ch <- mut.part.ch[,c(1,2,3,4,6,5)]
      mut.part.ch <- mut.part.ch[order(as.numeric(as.character(mut.part.ch$mutation))),]
      
      #exp.part.ch <- exp.part.d[which(exp.part$cohort==cohort),]
      #mut.part.ch <- mut.part.d[which(mut.part$cohort==cohort),]
      
      
      
      q <- expprint(exp.part.ch)
      #head(exp.part.ch)
      p <- mutprint.meta(mut.part.ch)
      p1 <- ggplot_gtable(ggplot_build(p))
      p2 <- ggplot_gtable(ggplot_build(q))
      maxWidth = unit.pmax(p1$widths[2:5], p2$widths[2:5])
      p1$widths[2:5] <- maxWidth
      p2$widths[2:5] <- maxWidth
      plot0 <- arrangeGrob(p2,p1, heights = c(2,3),
                           top=textGrob(paste(gene.name," inactivations",sep=""), gp=gpar(fontsize=14,font=8)))
      myplots[[i]] <<- plot0
      
    })
    i = i+1
  }
}else{
  myplot0 <- blankprint("inactivation and expression")
  inact_flag <- 1
}


# mutation cohort distribution
event <- "biallalic_deletion"
stat.matrix <- data.frame(t(gene.matrix))
stat.matrix$biallelic <-as.numeric(as.character(stat.matrix$biallelic))
distribution.matrix1 <- data.frame(aggregate(stat.matrix$biallelic,list(stat.matrix$cohort),sum))
colnames(distribution.matrix1) <- c("cohort",event)
distribution.matrix1 <- distribution.matrix1[distribution.matrix1$biallalic_deletion!=0,]
distribution.matrix1$cohort <- as.character(distribution.matrix1$cohort)
dh_no <- sum(distribution.matrix1$biallalic_deletion)

ba_entropy <- entropy(distribution.matrix1$biallalic_deletion)

if (dh_no != 0){
  r <- pieplot(distribution.matrix1,event,ba_entropy)
}else{
  r <- blankprint(paste(event," not detected",sep=""))
}


event <- "CNV_GAIN"
stat.matrix <- data.frame(t(gene.matrix))
stat.matrix$CNV_GAIN <-as.numeric(as.character(stat.matrix$CNV_GAIN))
distribution.matrix <- data.frame(aggregate(stat.matrix$CNV_GAIN,list(stat.matrix$cohort),sum))
colnames(distribution.matrix) <- c("cohort",event)
distribution.matrix <- distribution.matrix[distribution.matrix$CNV_GAIN!=0,]
distribution.matrix$cohort <- as.character(distribution.matrix$cohort)
gain_no <- sum(distribution.matrix$CNV_GAIN)

cng_entropy <- entropy(distribution.matrix$CNV_GAIN)

if (gain_no != 0){
s = pieplot(distribution.matrix,event,cng_entropy)

}else{
  s <-  blankprint(paste(event," not detected",sep=""))
}

# input file VI and VII
background_dh <- fread("/ibios/co02/chen/pcawg_2018/backgroud/backs2.txt",data.table = F)
background_gain <- fread("/ibios/co02/chen/pcawg_2018/backgroud/amp_overview_background.txt",data.table = F)
background_dh <- background_dh[,c(1,6)]
background_gain<- background_gain[,c(1,8)]

background_gain$toapmore <- as.numeric(as.character(background_gain$toapmore))
gain_rate = round(nrow(background_gain[background_gain$toapmore<gain_no,])/nrow(background_gain)*100,2)
j <- ggplot(background_gain,aes(x=toapmore))+geom_bar()+geom_vline(xintercept = gain_no,color="red")+
  ggtitle(paste(gene.name," has more copy number gain than ",gain_rate,"% genes",sep=""))+
  xlab("total CNV gain > 8")



background_dh$total_dh <- as.numeric(as.character(background_dh$total_dh))
gain_rate = round(nrow(background_dh[background_dh$total_dh<dh_no,])/nrow(background_dh)*100,2)
k <- ggplot(background_dh,aes(x=total_dh))+geom_bar()+geom_vline(xintercept = dh_no,color="red")+
  ggtitle(paste(gene.name," has more biallelic inactivations than ",gain_rate,"% genes",sep=""))+
  xlab("total biallelic inactivations")

pdf(paste(output.dir,"/",gene.name,".pdf",sep = ""),
    height = 24, width = 24)

grid.arrange(r,s,plot1,plot0,k,j,nrow=3,
             top = textGrob(paste(gene.name," functional mutations",sep=""), 
                       gp=gpar(fontsize=24,font=8))
             )
if (inact_flag == 0){
  pages <- floor(length(myplots)/4)
  for(i in 1:pages){
    k1 <- i*4-3
    plot(arrangeGrob(grobs=myplots[k1:(k1+3)],nrow=2))
  }
}
dev.off()


write.table(stat.matrix,file=paste(output.dir,"/",gene.name,".txt",sep = ""),
            quote = F, col.names = T, row.names = T, sep="\t")
