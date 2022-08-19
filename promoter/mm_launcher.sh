 
cohort_list="BLCA-US CESC-US ESAD-UK KIRP-US LIRI-JP OV-US PRAD-US UCEC-US BOCA-UK CLLE-ES GACA-CN LAML-KR LUAD-US PACA-AU READ-US BRCA-EU CMDI-UK GBM-US LAML-US LUSC-US PACA-CA SARC-US BRCA-UK COAD-US HNSC-US LGG-US MALY-DE PAEN-AU SKCM-US BRCA-US DLBC-US KICH-US LICA-FR ORCA-IN PBCA-DE STAD-US BTCA-SG EOPC-DE KIRC-US LIHC-US OV-AU PRAD-UK THCA-US LINC-JP PRAD-CA RECA-EU MELA-AU"

#cohort_list="PRAD-US SKCM-US BRCA-US"

for coh in ${cohort_list};do

job_mm=`qsub -j oe -o /home/hongc/PhD/promoters/log -v cohort=${coh} -N matrix_${coh} -l walltime=1:00:00  /home/hongc/PhD/promoters/code/st1_matrix_maker.sh| cut -d "." -f 1`

done
