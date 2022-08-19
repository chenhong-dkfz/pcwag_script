#cohs="SARC-US"
cohs="PBCA-DE PACA-AU PACA-CA PAEN-AU SARC-US OV-AU OV-US BRCA-EU BRCA-UK PRAD-UK PRAD-US EOPC-DE"
EO=/home/hongc/project/codes/enhancer_codes/output #output and error information


for coh in ${cohs};do

bed_path=/ibios/co02/chen/enhancer2/results
#bed_path=/ibios/co02/chen/promoter/results
bed_file=${bed_path}/${coh}.bed
snv_out=/ibios/co02/chen/enhancer2/results/analysis/${coh}_snv.bed
cnv_out=/ibios/co02/chen/enhancer2/results/analysis/${coh}_cnv.bed
sco_out=/ibios/co02/chen/enhancer2/results/analysis/${coh}_score.bed
#ref=/ibios/co02/chen/tmp/gencode.v19.genes.bed
ref=/ibios/co02/chen/enhancer2.bed
out_put=/ibios/co02/chen/enhancer2/results/analysis/${coh}_result.bed

rm ${out_put}
job_convert=`echo "python /home/hongc/project/codes/enhancer_codes/complete_code2/matrix_convert.py ${bed_file} ${snv_out} ${cnv_out}"|qsub -N ${coh} -j oe -o ${EO} -l walltime=00:15:00`
job_calc=`echo "Rscript /home/hongc/project/codes/enhancer_codes/complete_code2/scoring.r ${cnv_out} ${snv_out} ${sco_out}"|qsub -N ${coh}_scoring -j oe -o ${EO} -l walltime=5:00:00 -W depend="afterok:"${job_convert}`
#echo "Rscript /home/hongc/project/codes/enhancer_codes/scoring.r ${cnv_out} ${snv_out} ${sco_out}"|qsub -N ${coh}_scoring -j oe -o ${EO} -l walltime=00:30:00

#job_snv=`echo "sort -k4,4 -gr ${sco_out} | head -n 10 > ${sco_out}_snv"|qsub -N ${coh}_snv_trans -j oe -o ${EO} -l walltime=00:10:00 -W depend="afterok:"${job_calc}`
#job_cnv=`echo "sort -k5,5 -gr ${sco_out} | head -n 10 > ${sco_out}_cnv"|qsub -N ${coh}_cnv_trans -j oe -o ${EO} -l walltime=00:10:00 -W depend="afterok:"${job_calc}`
#job_mul=`echo "sort -k7,7 -gr ${sco_out} | head -n 10 > ${sco_out}_mul"|qsub -N ${coh}_mul_trans -j oe -o ${EO} -l walltime=00:10:00 -W depend="afterok:"${job_calc}`
#echo "bedtools intersect -a ${sco_out}_snv -b ${ref} -wao|head -n 10 |cut -f 1,2,3,4,12 >> ${out_put}"|qsub -N ${coh}_snv_map -j oe -o ${EO} -l walltime=00:10:00 -W depend="afterok:"${job_snv}
#echo "bedtools intersect -a ${sco_out}_cnv -b ${ref} -wao|head -n 10 |cut -f 1,2,3,5,12 >> ${out_put}"|qsub -N ${coh}_cnv_map -j oe -o ${EO} -l walltime=00:10:00 -W depend="afterok:"${job_cnv}
#echo "bedtools intersect -a ${sco_out}_mul -b ${ref} -wao|head -n 10 |cut -f 1,2,3,7,12 >> ${out_put}"|qsub -N ${coh}_mul_map -j oe -o ${EO} -l walltime=00:10:00 -W depend="afterok:"${job_mul}

done

