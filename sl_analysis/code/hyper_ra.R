library(data.table)
data <- fread("/abi/data/chen/project_august/mutation.csv",data.table=F)
data.new <- data[,-1]
data.new <- data.new[ , -which(names(data.new) %in% c("V444","V865"))]

# show high amplification, 0,1,3,4,5 -> 0, 2 -> 1
# show all amplifications, 0,3,4,5 -> 0, 1,2 -> 1
# show all mutations, 0,1,2 -> 0, 3,4,5 -> 1
# show all possible bas, 0,1,2,3 -> 0, 4,5 -> 1
# show all real bas, 0,1,2,3,4 -> 0, 5 -> 1


hyperdis <- function(a,b,mut_typ){
  
  if (mut_typ == "ham"){
    a <- replace(a,a==0|a==1|a==3|a==4|a==5,0)
    a <- replace(a,a==2,1)
    b <- replace(b,b==0|b==1|b==3|b==4|b==5,0)
    b <- replace(b,b==2,1)
  }
  
  else if (mut_typ == "aam"){
    a <- replace(a,a==0|a==3|a==4|a==5,0)
    a <- replace(a,a==2,1)
    b <- replace(b,b==0|b==3|b==4|b==5,0)
    b <- replace(b,b==2,1)
  }
  
  else if (mut_typ == "mut"){
    a <- replace(a,a==0|a==1|a==2,0)
    a <- replace(a,a==3|a==4|a==5,1)
    b <- replace(b,b==0|b==1|b==2,0)
    b <- replace(b,b==3|b==4|b==5,1)
  }
  
  else if (mut_typ == "pba"){
    a <- replace(a,a==0|a==1|a==2|a==3,0)
    a <- replace(a,a==4|a==5,1)
    b <- replace(b,b==0|b==1|b==2|b==3,0)
    b <- replace(b,b==4|b==5,1)
  }
  
  else if (mut_typ == "rba"){
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

testdf <- data.new
#data[data$gene_name=="PTEN",1:5]
target.gene <- data.new[24605,] #
#mut0 <- apply(testdf, 1, hyperdis,b=target.gene,mut_typ="mut")
#save(mut0,file="/abi/data/chen/mut1.RData")
#pba0 <- apply(testdf, 1, hyperdis,b=target.gene,mut_typ="pba")
#save(mut0,file="/abi/data/chen/pba1.RData")
rba0 <- apply(testdf, 1, hyperdis,b=target.gene,mut_typ="rba")
save(rba0,file="/abi/data/chen/rba1.RData")