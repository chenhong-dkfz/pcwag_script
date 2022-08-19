 
cohs="PBCA-DE PACA-AU PACA-CA PAEN-AU SARC-US OV-AU OV-US BRCA-EU BRCA-UK PRAD-UK PRAD-US EOPC-DE" #applied cohorts
#cohs="BRCA-EU BRCA-UK"
code_dir=/home/hongc/CaSINo/code
analysis_scripts_dir=${code_dir}/analysis_scripts
output_dir=/home/hongc/test/log


for coh in ${cohs};do

mkdir -p /ibios/co02/chen/new_casino2/matrix3/enhancer_matrix/${coh}

#qsub -N ${coh}_parser -v coh=${coh} -j oe -o ${output_dir} ${analysis_scripts_dir}/vcf2bed.sh -l walltime=2:00:00

#bed_file=/ibios/co02/chen/new_casino2/results/enhancer/${coh}.bed
#snv_out=/ibios/co02/chen/new_casino2/matrix/enhancer_matrix/${coh}/${coh}_snv.bed
#cnv_out=/ibios/co02/chen/new_casino2/matrix/enhancer_matrix/${coh}/${coh}_cnv.bed

#echo "python ${analysis_scripts_dir}/matrix_convert.py yes ${bed_file} ${snv_out} ${cnv_out}" | qsub -N ${coh}_cvt -v coh=${coh} -j oe -o ${output_dir} -l walltime=00:15:00

for turn in 1 2 3 4 5 6 7 8 9 10;do

#turn='1'
#bed_file=/ibios/co02/chen/new_casino3/results/promoter_randoms/${coh}.bed
bed_file=/ibios/co02/chen/new_casino3/results/promoter_randoms/${coh}_${turn}.bed
snv_out=/ibios/co02/chen/new_casino3/matrix/promoter_matrix/${coh}/${coh}_t${turn}_snv.bed
cnv_out=/ibios/co02/chen/new_casino3/matrix/promoter_matrix/${coh}/${coh}_cnv.bed

mkdir -p /ibios/co02/chen/new_casino3/matrix/promoter_matrix/${coh}

#job_convert=`qsub -N ${coh}_parser -v coh=${coh},turn=${turn} -j oe -o ${output_dir} ${analysis_scripts_dir}/vcf2bed_ran.sh -l walltime=2:00:00`

#echo "python ${analysis_scripts_dir}/matrix_convert.py no ${bed_file} ${snv_out} ${cnv_out}" | qsub -N ${coh}_${turn}_cvt -v coh=${coh} -j oe -o ${output_dir} -l walltime=00:15:00 -W depend="afterok:"${job_convert} 
echo "python ${analysis_scripts_dir}/matrix_convert.py no ${bed_file} ${snv_out} ${cnv_out}" | qsub -N ${coh}_${turn}_cvt -v coh=${coh} -j oe -o ${output_dir} -l walltime=00:15:00 

done


done