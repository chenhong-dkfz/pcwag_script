library(data.table)
library(ComplexHeatmap)
library(dplyr)
library(easyGgplot2)
library(grid)
library(gridExtra)
library(ggplot2)

data <- fread("/abi/data/chen/project_august/mutation2.csv",data.table=F)


blankprint <- function(gn,title){
  bplot <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank() + 
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = 14,hjust = 0.5))+
    ggtitle(paste(gn," ",title,sep=""))
  
  return(bplot)
}




make_plot <- function(gene1,gene2,gene3){

if(gene3 != "none"){data.new <- data[data$gene_name==gene1|data$gene_name==gene2|data$gene_name==gene3,]}else{
  data.new <- data[data$gene_name==gene1|data$gene_name==gene2,]
}
rownames(data.new) <- data.new$gene_name
data.new <- data.new[,-1]
#data.new <- data.new[ , -which(names(data.new) %in% c("V444","V865"))]

data.new <- replace(data.new,data.new==5,"BI")
data.new <- replace(data.new,data.new==4,"BI")
data.new <- replace(data.new,data.new==3,"MUT")
data.new <- replace(data.new,data.new==2,"AMP")
data.new <- replace(data.new,data.new==1,"AMP")
data.new <- replace(data.new,data.new==0,"")

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  BI = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue", col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red", col = NA))
  },
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#008000", col = NA))
  }
)

col = c("MUT" = "#008000", "AMP" = "red", "BI" = "blue")

mat <- data.new
oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "OncoPrint for ICGC pan-cancer, codeletion exclusive genes",
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "BI", "MUT"), 
                                      labels = c("Amplification", "Deep deletion", "Mutation")))


cohort.data <- fread("/abi/data/chen/project_august/pcawg_ta_project.txt",data.table = F)
t_mat <- data.frame(t(mat))
t_mat$pid <- rownames(t_mat)
t_mat_s <- merge(t_mat,cohort.data,by.x="pid",by.y="tumor_aliquot_id_wgs")
mylist <- split(t_mat_s, t_mat_s$project_code)


myplots <- list()
mats <- list()
for(i in 1:length(mylist)){
  print(i)
  t_mat_i <- mylist[[i]]
  cohort.id <- t_mat_i$project_code[1]
  t_mat_i <- t_mat_i[,-which(names(t_mat_i) %in% c("project_code"))]
  mat_i <- t(t_mat_i)
  colnames(mat_i) <- mat_i[1,]
  mats[[i]] <- mat_i[-1,]
  
  if (length(grep("AMP", mats[[i]]))==0 & length(grep("BI", mats[[i]]))==0 & length(grep("MUT", mats[[i]]))==0 ){
    myplots[[i]] <- blankprint(cohort.id," has no variants")
  }else{
  
  ht <- oncoPrint(mats[[i]], get_type = function(x) strsplit(x, ";")[[1]],
            alter_fun = alter_fun, col = col, 
            column_title = paste0(cohort.id," , codeletion exclusive genes"),
            heatmap_legend_param = list(title = "Alternations", at = c("AMP", "BI", "MUT"), 
                                        labels = c("Amplification", "Deep deletion", "Mutation")))
  myplots[[i]] <- grid.grabExpr({ draw(ht)})
  }
}


pdf(paste0("/home/hongc/test/mutual_exclusive/",gene1,"_",gene2,".pdf"),width = 20)
#ggplot2.multiplot(plotlist = myplots, cols = 3)
pages <- floor(length(myplots)/4)
for(i in 1:pages){
  k1 <- i*4-3
  grid.newpage()
  x<-arrangeGrob(grobs=myplots[k1:(k1+3)],nrow=2)
  grid.draw(x)
}
dev.off()
}

make_plot("PTEN","SMAD4","none")
make_plot("ARID1A","ARID1B","none")
make_plot("PTEN","PLK4")
make_plot("PTEN","BRCA1","BRCA2")
make_plot("VHL","FMN1","none")

