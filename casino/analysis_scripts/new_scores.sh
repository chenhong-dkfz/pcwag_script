cohs="PBCA-DE PACA-AU PAEN-AU SARC-US OV-AU OV-US BRCA-EU BRCA-UK PRAD-UK PRAD-US EOPC-DE" #applied cohorts
#cohs="PACA-CA"
df0=/ibios/co02/chen/new_casino3/matrix/promoter_matrix #random path
df=/ibios/co02/chen/new_casino2/matrix/promoter_matrix #real path
output_dir=/home/hongc/CaSINo/code/job_output
analysis_scripts_dir=/home/hongc/CaSINo/code/analysis_scripts

for coh in ${cohs};do

snv_out=${df}/${coh}/${coh}_snv.bed
cnv_out=${df}/${coh}/${coh}_cnv.bed
sco=/ibios/co02/chen/new_casino3/reports/promoter_${coh}.txt
job_calc=`echo "Rscript ${analysis_scripts_dir}/scores_june.r ${snv_out} ${cnv_out} ${df0} ${sco} ${coh}"|qsub -N ${coh}_scoring -j oe -o ${output_dir} -l walltime=01:40:00`

done
