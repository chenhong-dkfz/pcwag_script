library(data.table)
library(ggplot2)
library(stringr)

# snv and indel, 2961 samples, 378 empty files, 2574 samples
# cnv, 2778 samples, 43 empty files, 2735 samples
# sv, 2883 samples, 218 empty files, 2665 samples
# rna-seq, 1151 samples


setwd("/abi/data/chen/project_august/core")
fhis <- fread("pcawg_specimen_histology_August2016_v9.tsv",data.table = F)
fcli <- fread("pcawg_donor_clinical_August2016_v9.tsv",data.table = F)

f1 <- fcli[,c(1,2,6,11,12,18)]
colnames(f1)[1] <- "donor_unique_id"
fhis <- fhis[fhis$specimen_library_strategy!="RNA-Seq",]
f2 <- fhis[,c(8,12,17,20)] # tier 4 has different values in same sample
f3 <- f2[!duplicated(f2), ]
fall <- merge(f1,f3,by="donor_unique_id")
fall <- fall[fall$project_code!="",]
idsplit <- data.frame(fall,do.call(rbind,str_split(fall$donor_unique_id,"::")))
fall$donor_unique_id <- idsplit$X2

fmap <- fread("PCAWG_id_table_donor_tumor_wgs.tsv",data.table = F)
fall2 <- merge(fall,fmap,by.x = "donor_unique_id",by.y = "PID")




fa <- fread("res.txt",data.table = F)
fb <- fread("re.txt",data.table = F,header = F)
fc <- data.frame(fa,do.call(rbind,str_split(fa$V2,"/")))
fd <- data.frame(fc,do.call(rbind,str_split(fc$X7,"\\.")))
fe <- fd[,c(1,8,10)]
names(fe) <- c("num","type","id")
ff <- merge(fe,fb,by.x="id",by.y="V3")
fg <- ff[ff$num!=0 & ff$type=="snv_mnv_del",]
fh <- ff[,c(3,5)]
fi <- fg[,c(4,5,1)]
fi <- unique(fi)

# remove duplicated from fmap
# add fixed fi, filtered by snv

fi$V2 <- substr(fi$V2,1,nchar(fi$V2)-2)
fj <- fmap[!fmap$PID %in% fb$V2,]
colnames(fi) <- colnames(fj)
fk <- rbind(fj,fi)

fall2 <- merge(fall,fk,by.x = "donor_unique_id",by.y = "PID")
#fall[!fall$donor_unique_id %in% fall2$donor_unique_id,]
#fall2[!fall2$donor_unique_id %in% fk$PID,]
#f0 <- fall[duplicated(fall$donor_unique_id), ]
write.table(fall2,file="sample_info_tumor_aliquot.txt",col.names = TRUE, row.names = F,quote = F,sep="\t")


stats <- fread("sample_stats.txt",data.table = F)
snv <- fread("sumsnv.txt",data.table = F)
snv$V2 <- substr(snv$V2,1,nchar(snv$V2)-4)
# q1 <- snv[snv$V1!=0,]
# q2 <- stats[stats$V2!=0,]

fall3 <- cbind(fall2,stats)
fall3 <- fall3[fall3$tumor_aliquot_id_wgs == fall3$V1,]
fall3 <- fall3[,c(1,2,3,4,5,6,7,8,9,10,11,13,14,15)]
colnames(fall3)[12] <- "SNV_num"
colnames(fall3)[13] <- "CNV_num"
colnames(fall3)[14] <- "SV_num"
fall3$SNV <- 1
fall3[fall3$SNV_num==0,15] <- 0
fall3$CNV <- 1
fall3[fall3$CNV_num==0,16] <- 0
fall3$SV <- 1
fall3[fall3$SV_num==0,17] <- 0

expedstar <- fread("exped_star.txt",data.table = F,header = F)
rnamap <- fread("ta_star2.txt",data.table = F)
rnamap <- rnamap[rnamap$tumor_rna_seq_star_alignment_bam_file_name %in% expedstar$V1,]
fall4 <- merge(fall3,rnamap,by.x = "tumor_aliquot_id_wgs",by.y = "Tumor_wgs_aliquot_id",all.x = T)
fall4[duplicated(fall4$tumor_aliquot_id_wgs), ]
fall4$RNA <- 1
fall4[is.na(fall4$tumor_rna_seq_star_alignment_bam_file_name)== TRUE,19] <- 0
write.table(fall5,file="pcawg_sample_info_v2.txt",col.names = TRUE, row.names = F,quote = F,sep="\t")

## summary

table(fall4$project_code)
table(fall4$organ_system)
table(fall4$tumour_stage)
table(fall4$donor_sex)
table(fall4$donor_age_at_diagnosis)
table(fall4$donor_wgs_included_excluded)

table(fall5$SNV)
table(fall5$CNV)
table(fall5$SV)
table(fall5$RNA)


## visualization
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
# SNV
fall5 <- fall4[fall4$donor_wgs_included_excluded == "Included",]
fall5$SNV <- as.factor(fall5$SNV)
snv_proj <- ggplot(fall5[order(fall5$SNV),],aes(y=1,x=project_code,fill=SNV))+geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle("samples with SNV")
snv_org <- ggplot(fall5[order(fall5$SNV),],aes(y=1,x=organ_system,fill=SNV))+geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle("samples with SNV")

# CNV
#fall5 <- fall4
fall5$CNV <- as.factor(fall5$CNV)
cnv_proj <- ggplot(fall5[order(fall5$CNV),],aes(y=1,x=project_code,fill=CNV))+geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle("samples with CNV")
cnv_org <- ggplot(fall5[order(fall5$CNV),],aes(y=1,x=organ_system,fill=CNV))+geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle("samples with CNV")

# CNV
#fall5 <- fall4
fall5$SV <- as.factor(fall5$SV)
sv_proj <- ggplot(fall5[order(fall5$SV),],aes(y=1,x=project_code,fill=SV))+geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle("samples with SV")
sv_org <- ggplot(fall5[order(fall5$SV),],aes(y=1,x=organ_system,fill=SV))+geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle("samples with SV")

# CNV
#fall5 <- fall4
fall5$RNA <- as.factor(fall5$RNA)
rna_proj <- ggplot(fall5[order(fall5$RNA),],aes(y=1,x=project_code,fill=RNA))+geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle("samples with RNA")
rna_org <- ggplot(fall5[order(fall5$RNA),],aes(y=1,x=organ_system,fill=RNA))+geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle("samples with RNA")

glist1 <- grid.arrange(snv_proj,snv_org,cnv_proj,cnv_org,sv_proj,sv_org,rna_proj,rna_org)

pdf("samples_info.pdf")
plot(glist1)
dev.off()






