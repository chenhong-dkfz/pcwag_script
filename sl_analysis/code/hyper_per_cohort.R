library(data.table)
data <- fread("/abi/data/chen/project_august/mutation.csv",data.table=F)
data.new <- data[,-1]
data.new <- data.new[ , -which(names(data.new) %in% c("V444","V865"))]


hyperdis <- function(a,b,mut_typ){
  
  if (mut_typ == "ham"){
    a <- replace(a,a==0|a==1|a==3|a==4|a==5,0)
    a <- replace(a,a==2,1)
    b <- replace(b,b==0|b==1|b==3|b==4|b==5,0)
    b <- replace(b,b==2,1)
  }else if (mut_typ == "aam"){
    a <- replace(a,a==0|a==3|a==4|a==5,0)
    a <- replace(a,a==2,1)
    b <- replace(b,b==0|b==3|b==4|b==5,0)
    b <- replace(b,b==2,1)
  }else if (mut_typ == "mut"){
    a <- replace(a,a==0|a==1|a==2,0)
    a <- replace(a,a==3|a==4|a==5,1)
    b <- replace(b,b==0|b==1|b==2,0)
    b <- replace(b,b==3|b==4|b==5,1)
  }else if (mut_typ == "pba"){
    a <- replace(a,a==0|a==1|a==2|a==3,0)
    a <- replace(a,a==4|a==5,1)
    b <- replace(b,b==0|b==1|b==2|b==3,0)
    b <- replace(b,b==4|b==5,1)
  }else if (mut_typ == "rba"){
    a <- replace(a,a==0|a==1|a==2|a==3|a==4,0)
    a <- replace(a,a==5,1)
    b <- replace(b,b==0|b==1|b==2|b==3|b==4,0)
    b <- replace(b,b==5,1)
  }
  
  mut1 <- sum(a==1)
  mut2 <- sum(b==1)
  inter <- sum((a+b)==2)
  len <- length(a)
  #print(mut1)
  #print(mut2)
  p.value <- phyper(inter,mut1,len-mut1,mut2,lower.tail = T)
}

testdf <- data.new[1:10,]
#data[data$gene_name=="PTEN",1:5]
target.gene <- data.new[24605,] #

cohort.data <- fread("/abi/data/chen/project_august/pcawg_ta_project.txt",data.table = F)

t_mat <- data.frame(t(data.new))
#t_mat <- data.frame(t(testdf))


t_mat$pid <- rownames(t_mat)
t_mat_s <- merge(t_mat,cohort.data,by.x="pid",by.y="tumor_aliquot_id_wgs")
mylist <- split(t_mat_s, t_mat_s$project_code)
pvalue <- list()
cohorts <- list()

for (i in 1:length(mylist)){
  print(i)
  cohort.id <- mylist[[i]]$project_code[1]
  cohorts[[i]] <- cohort.id
  rownames(mylist[[i]]) <- mylist[[i]]$pid
  mylist[[i]] <- mylist[[i]][,-which(names(mylist[[i]]) %in% c("project_code"))]
  mat_i <- data.frame(t(mylist[[i]]),check.names = F)
  mat_i <- mat_i[-1,]
  mat_i[] <- lapply(mat_i, function(x) as.numeric(as.character(x)))
  #target.g <- mat_i[24605,]
  target.g <- mat_i[2,]
  
  pvalue[[i]] <- apply(mat_i, 1, hyperdis,b=target.g,mut_typ="mut")
  
  }
save(pvalue,file="/abi/data/chen/mut_per_cohort.RData")
save(cohorts,file="/abi/data/chen/cohorts.RData")

