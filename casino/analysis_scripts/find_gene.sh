cohs="PBCA-DE PACA-AU PACA-CA PAEN-AU SARC-US OV-AU OV-US BRCA-EU BRCA-UK PRAD-UK PRAD-US EOPC-DE" #applied cohorts
result_path=/ibios/co02/chen/casino/promoter/results


for coh in ${cohs};do
  awk '/\./' ${result_path}/${coh}.bed > ${result_path}/${coh}_r.bed
  #awk 'NR != 14466' ${result_path}/${coh}.bed > ${result_path}/${coh}_r.bed
  bedtools intersect -a ${result_path}/${coh}_r.bed -b ${result_path}/analysis/selected_candidates_${coh}.txt -wao > ${result_path}/${coh}_snv.bed
  

#python /home/hongc/CaSINo/code/analysis_tools/find_gene.py /ibios/co02/chen/casino/promoter/results/analysis/selected_candidates_${coh}.txt
done