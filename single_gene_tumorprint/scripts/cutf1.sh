cohort_list="BLCA-US CESC-US ESAD-UK KIRP-US LIRI-JP OV-US PRAD-US UCEC-US BOCA-UK CLLE-ES GACA-CN LAML-KR LUAD-US PACA-AU READ-US BRCA-EU CMDI-UK GBM-US LAML-US LUSC-US PACA-CA SARC-US BRCA-UK COAD-US HNSC-US LGG-US MALY-DE PAEN-AU SKCM-US BRCA-US DLBC-US KICH-US LICA-FR ORCA-IN PBCA-DE STAD-US BTCA-SG EOPC-DE KIRC-US LIHC-US OV-AU PRAD-UK THCA-US LINC-JP PRAD-CA RECA-EU MELA-AU PAEN-IT"

for cohort in ${cohort_list};do
    cut -f 1 /ibios/co02/chen/pcawg_2018/pid_list/whitelist/${cohort}.txt > /ibios/co02/chen/pcawg_2018/pid_list/whitelist0/${cohort}.txt
done