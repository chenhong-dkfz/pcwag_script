#install.packages("writexl")
#library(writexl)
library(data.table)
options(java.parameters = "-Xmx1024m")
library(XLConnect)

#library(XLConnect)

gene = "ALC1"

p.table = fread(paste0("/abi/data/chen/project_august/mutation25/cohortwise/full2/",tolower(gene),"_hyper.txt"),
                  data.table = F)
o.table = fread(paste0("/abi/data/chen/project_august/mutation25/cohortwise/full2/",tolower(gene),"_oddsratio.txt"),
                data.table = F)
cgc.table <- fread("/home/hongc/PhD/basic/cgc.txt",data.table = F)

mutss <- c("mut_mut","mut_nam","mut_pam","mut_pa","mut_ra",
           "nam_mut","nam_nam","nam_pam","nam_pa","nam_ra",
           "pam_mut","pam_nam","pam_pam","pam_pa","pam_ra",
           "ra_mut","ra_nam","ra_pam","ra_pa","ra_ra",
           "pa_mut","pa_nam","pa_pam","pa_pa","pa_ra")

p.slide <- p.table[,c(1,2)]
p.slide$cgc <- F
p.slide[p.slide$gene_name %in% cgc.table$`Gene Symbol`,3] <- T

wb <- loadWorkbook(paste0("/home/hongc/test/",gene,".xls"),
 create = TRUE)

for (mut in mutss){
  if(gene=="ALC1"){
    gene_alias <- "CHD1L"
  }
  stats <- fread(paste0("/abi/data/chen/project_august/mutation25/cohortwise/",gene_alias,"/",mut,".txt"),
                        data.table=F)
  
  colnames(stats) <- c("sample_size",
                       paste0(gene,"_mutations"),
                       "gene_mutation",
                       "overlaps")
  
  p.sheet <- cbind(p.slide,p.table[,mut],o.table[,mut])
  colnames(p.sheet)[4] <-  "p"
  colnames(p.sheet)[5] <-  "OR"
  p.sheet <- cbind(p.sheet,stats)
  p.sheet[p.sheet$OR=="Inf",5] <- 1
  #write.xlsx2(p.sheet, file = paste0("/home/hongc/test/",gene,".xls"), sheetName=mut, append=TRUE)
  
  #createSheet(wb, name = mut)
  #writeWorksheet(wb, p.sheet, sheet = mut)
  #saveWorkbook(wb)
  
  write.table(p.sheet,
              file=paste0("/abi/data/chen/project_august/mutation25/cohortwise/result/",gene,"/",mut,".txt"),
              col.names = T, row.names = F, quote = F, sep = "\t")
}
