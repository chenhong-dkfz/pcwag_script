library("cli")
library("writexl")
library(data.table)

setwd("/home/hongc/test/result_subcohort2")


muts <- c("mut","nam","pam","pa","ra")
mutss <- NULL
i <- 1
for (mut in muts){
  for (mut2 in muts){
    st <- paste0(mut,"_",mut2)
    mutss[i] <- st
    i <- i + 1
  }
}

for(mut in mutss){
  setwd(paste0("/home/hongc/test/result_subcohort2/",mut))
  
  FileNameExport <- paste0("/home/hongc/test/result_subcohort2/",mut)
  
  temp = list.files(pattern="*.txt")
  myfiles = lapply(temp, fread,data.table=F)
  for(i in 1:length(myfiles)){
   myfiles[[i]] <- myfiles[[i]][,c(1,4,7,8,9,10,11,12)]
  }
  
  legend_descriptor <- NULL
  legend_descriptor <- data.frame(
    label=c(
      ## selected subcohort
      "dataset",
      "histology",
      "gene_of_interest",
      ## the names of the tables (xlsx sheet names corresponding to the 25 comparisons), e.g.:
      "mut_hiCo",
      "mut_medCo",
      "mut_lowCo",
      "amp_hi",
      "amp_norm",
      ## table fields
      "chr",
      "gene_name",
      "cgc",
      "p",
      "OR",
      "sample_size",
      #         "PTEN_mutations",
      "CHD1L_mutations",
      "gene_mutation",
      "overlaps"
    ),
    description=c(
      ## selected subcohort
      "BRCA-EU, BRCA-US, ...",
      "breast, breast, ...",
      #         "PTEN",
      "CHD1L",
      ## the names of the tables (xlsx sheet names corresponding to the 25 comparisons), e.g.:
      "bi-allelic mutations only",
      "bi-allelic + potential bi-allelic mutations",
      "bi-allelic + potential bi-allelic + mono-allelic mutations + ",
      "copy number gains >8",
      "copy number gains >4",
      ## table fields
      "chromosome the gene is annotated to",
      "official gene symbol by HUGO nomenclature",
      "cancer gene census annotation",
      "p-value of hypergeometric test",
      "effect strength (odds ratio)",
      "total number of samples in the selected subcohort(s)",
      #         "count of samples with PTEN mutation",
      "count of samples with CHD1L mutation",
      "count of samples with mutation in other gene",
      "count of samples with mutation in both genes"
    ),
    stringsAsFactors=FALSE
  )
  
  names_vector  <- gsub(paste0("_",mut,".txt"),"",temp)
  datalist <- myfiles
  stopifnot(length(datalist) == length(names_vector))
  names(datalist) <- names_vector
  #datalist[["legend"]] <- legend_descriptor
  ## re-order legend to beginning
  #datalist <- datalist[c("legend", names_vector)]
  cli::cat_line("Writing ", length(datalist), " tables to file 'PTEN_", FileNameExport, ".txt'.")
  write_xlsx(x=datalist, path=paste0("/home/hongc/test/result_101219/PTEN/PTEN_",mut,".xlsx"), col_names=TRUE)
  
}
