cohs="PBCA-DE PACA-AU PACA-CA PAEN-AU SARC-US OV-AU OV-US BRCA-EU BRCA-UK PRAD-UK PRAD-US EOPC-DE"
path=/home/hongc/CaSINo/code/small_tools/sv_test
path2=/ibios/co02/chen/casino_sv
for coh in ${cohs};do
job1=`echo "python ${path}/matrix_convert.py ${path}/result/enhancer/${coh}/${coh}.bed ${path2}/temp/${coh}.bed"|qsub -j oe -o ${path}/log -N mc_${coh} -l walltime=00:20:00`
job2=`echo "Rscript ${path}/scoring_sv.r ${path2}/temp/${coh}.bed ${path2}/result/enhancer/${coh}_scored.bed"|qsub -j oe -o ${path}/log -N sc_${coh} -l walltime=00:20:00 -W depend="afterok:"${job1}`
done
