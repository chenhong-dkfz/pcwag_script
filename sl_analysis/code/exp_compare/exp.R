require(data.table)
require(reshape2)

fin <- fread("/abi/data/chen/project_august/exp_comp/im.txt",data.table = F)
fin.t <- data.frame(t(fin))
gene1 <- as.character(fin.t[1,1])
gene2 <- as.character(fin.t[1,2])
fin.t$pid <- rownames(fin.t)
fin.t <- fin.t[-1,]
colnames(fin.t) <- c(gene1,gene2,"pid")
cohort.data <- fread("/abi/data/chen/project_august/pcawg_ta_project.txt",data.table = F)
pid.conv <- fread("/home/hongc/PhD/basic/tumor_rna_new_uniq.txt",data.table = F)
fin.t <- merge(fin.t,pid.conv,by.x="pid",by.y = "tumor_rna_seq_star_alignment_bam_file_name")
fin.t <- merge(fin.t,cohort.data,by.x="Tumor_wgs_aliquot_id",by.y="tumor_aliquot_id_wgs")
fin.t$CHD1 <- as.numeric(as.character(fin.t$CHD1))
fin.t$PTEN <- as.numeric(as.character(fin.t$PTEN))
coh.his <- fread("/home/hongc/PhD/basic/cohort_histology.tsv",data.table = F,header = F)
colnames(coh.his) <- c("project_code","histology")
fin.t <- merge(fin.t,coh.his,by="project_code")


# plot 1 
ggplot(fin.t,aes(y=CHD1,x=PTEN,color=histology))+geom_point()+facet_wrap(~histology)

# plot 2
fin.t <- fin.t[order(fin.t$PTEN),]
fin.t1 <- fin.t
fin.m1 <- melt(fin.t1)
fin.m1 <- fin.m1[order(fin.m1$project_code),]
ggplot(fin.m1,aes(y=value,x=pid,colour=variable)) +
  geom_point() +
  facet_wrap(~project_code) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggplot(fin.t1[fin.t1$project_code=="PRAD-US",],aes(y=CHD1,x=PTEN))+geom_point()

