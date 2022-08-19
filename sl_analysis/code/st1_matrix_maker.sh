#parameter setting
source_dir="/ibios/co02/chen/pcawg_2018"
code_dir="/home/hongc/PhD/TumorPrint/single_gene/code"
roi=/ibios/co02/chen/pcawg_2018/regions/coding.txt
#cohort_list="BLCA-US CESC-US ESAD-UK KIRP-US LIRI-JP OV-US PRAD-US UCEC-US BOCA-UK CLLE-ES GACA-CN LAML-KR LUAD-US PACA-AU READ-US BRCA-EU CMDI-UK GBM-US LAML-US LUSC-US PACA-CA SARC-US BRCA-UK COAD-US HNSC-US LGG-US MALY-DE PAEN-AU SKCM-US BRCA-US DLBC-US KICH-US LICA-FR ORCA-IN PBCA-DE STAD-US BTCA-SG EOPC-DE KIRC-US LIHC-US OV-AU PRAD-UK THCA-US LINC-JP PRAD-CA RECA-EU MELA-AU"

#cohort_list="PRAD-UK"

#information parser(snv,sv,cnv)

# python snv_parser.py
# python sv_parser.py
# python cnv_parser.py
# python region_of_interests_preparation.py


#matrix per gene

module load bedtools/2.16.2

mkdir -p ${source_dir}/result/matrix/coding

#for cohort in ${cohort_list};do
    pid_file="/ibios/co02/chen/pcawg_2018/pid_list/whitelist/${cohort}.txt"
    while read line; do
        mkdir -p ${source_dir}/result/${cohort}/${line}
        
        #echo ${line}
        
        bedtools intersect -a ${roi} -b ${source_dir}/snv_exon2/${line}.txt -wao > ${source_dir}/result/${cohort}/${line}/snv_tmp.bed
        python ${code_dir}/merge_same_positions.py ${source_dir}/result/${cohort}/${line}/snv_tmp.bed ${source_dir}/result/${cohort}/${line}/snv.bed 
        cut -f 5 ${source_dir}/result/${cohort}/${line}/snv.bed > ${source_dir}/result/${cohort}/${line}/snv_tmp.bed
        mv ${source_dir}/result/${cohort}/${line}/snv_tmp.bed ${source_dir}/result/${cohort}/${line}/snv.bed
        
        bedtools intersect -a ${roi} -b ${source_dir}/cnv/${line}.txt -wao > ${source_dir}/result/${cohort}/${line}/cnv_tmp.bed
        python ${code_dir}/merge_same_positions.py ${source_dir}/result/${cohort}/${line}/cnv_tmp.bed ${source_dir}/result/${cohort}/${line}/cnv.bed 
        cut -f 5 ${source_dir}/result/${cohort}/${line}/cnv.bed > ${source_dir}/result/${cohort}/${line}/cnv_tmp.bed
        mv ${source_dir}/result/${cohort}/${line}/cnv_tmp.bed ${source_dir}/result/${cohort}/${line}/cnv.bed
        
        bedtools intersect -a ${roi} -b ${source_dir}/sv/${line}.txt -wao > ${source_dir}/result/${cohort}/${line}/sv_tmp.bed
        python ${code_dir}/merge_same_positions.py ${source_dir}/result/${cohort}/${line}/sv_tmp.bed ${source_dir}/result/${cohort}/${line}/sv.bed 
        cut -f 5 ${source_dir}/result/${cohort}/${line}/sv.bed > ${source_dir}/result/${cohort}/${line}/sv_tmp.bed
        mv ${source_dir}/result/${cohort}/${line}/sv_tmp.bed ${source_dir}/result/${cohort}/${line}/sv.bed
    
        
    paste -d ';' ${source_dir}/result/${cohort}/${line}/snv.bed ${source_dir}/result/${cohort}/${line}/cnv.bed ${source_dir}/result/${cohort}/${line}/sv.bed > ${source_dir}/result/${cohort}/${line}/roi.bed
    sed -i -- 's/\(;;\|;;;\)/;/g' ${source_dir}/result/${cohort}/${line}/roi.bed
    sed -i -- 's/\(^;\|;$\)//g' ${source_dir}/result/${cohort}/${line}/roi.bed
    sed -i -- 's/\(\.;\|;\.\)//g' ${source_dir}/result/${cohort}/${line}/roi.bed
    sed -i '1s/^/'${line}'\n/' ${source_dir}/result/${cohort}/${line}/roi.bed 
   
    
    done<${pid_file}
    paste ${roi} `ls  ${source_dir}/result/${cohort}/*/roi.bed` > ${source_dir}/result/matrix/coding/${cohort}2.bed
    sort -k4,4 ${source_dir}/result/matrix/coding/${cohort}2.bed > ${source_dir}/result/matrix/coding/${cohort}.bed


#oncoprint
 
#['nonsynonymous SNV', 'synonymous SNV', 'stopgain', 'nonframeshift substitution', 'frameshift insertion', 'frameshift deletion', 'stoploss', 'unknown', 'nonframeshift deletion', 'nonframeshift insertion', 'frameshift substitution']
