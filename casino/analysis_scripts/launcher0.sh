 
cohs="PBCA-DE PACA-AU PACA-CA PAEN-AU SARC-US OV-AU OV-US BRCA-EU BRCA-UK PRAD-UK PRAD-US EOPC-DE" #applied cohorts
#cohs="PBCA-DE PACA-AU PAEN-AU SARC-US OV-AU OV-US BRCA-EU BRCA-UK PRAD-UK PRAD-US EOPC-DE"
#cohs="PACA-CA"
code_dir=/home/hongc/CaSINo/code
analysis_scripts_dir=${code_dir}/analysis_scripts
output_dir=/home/hongc/test/log


for coh in ${cohs};do

qsub -v coh=${coh} -l walltime=02:00:00 -N ${coh}_parsing -j oe -o ${output_dir} ${analysis_scripts_dir}/vcf2bed_snv.sh

done 
