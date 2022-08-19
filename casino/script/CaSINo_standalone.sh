#parameter

#cohorts
#cohs="PBCA-DE PACA-AU PACA-CA PAEN-AU SARC-US OV-AU OV-US BRCA-EU BRCA-UK PRAD-UK PRAD-US EOPC-DE" #applied cohorts
#cohs="PRAD-UK" #cohort for test only
#cohs=`cat /ibios/co02/chen/data_march_2016/cohort_list.txt`
cohs="BLCA-US SKCM-US THCA-US"

#directory which contains pid lists per cohort
pid_path=/home/hongc/CaSINo/coh_pid 
#pid_path=/ibios/co02/chen/data_march_2016/variants_per_cohorts/temp

#roi file is a file with 4 columns,chrom,start,end,information and ends with column names
#roi_file=${pid_path}/enhancer.txt  #enhancer
roi_file=${pid_path}/promoter.txt  #promoter
#roi_file=${pid_path}/coding.bed #coding part

#main directory of codes
code_dir=/home/hongc/CaSINo/code
analysis_scripts_dir=${code_dir}/analysis_scripts

#directory for results
result_path=/ibios/co02/chen/work/08_07_2016

#result directory for random real output
result_dir=${result_path}/enhancer/results
#result_dir=${result_path}/promoter/results
#result_dir=/ibios/co02/chen/casino/results

#directory which contains intermediate files per cohort
#mut_path=/ibios/co02/chen/casino/intermediate
#mut_path=${result_path}/promoter/intermediate
mut_path=${result_path}/enhancer/intermediate
#job reports
output_dir=/home/hongc/CaSINo/code/job_output



#IMPORTANT!#
#you have to modify plotqq.r for the parsed RANDOMIZED MUTATION DATA because I cannot find original MUTATION FILE any more...#

######
#jobs#
######

mkdir -p ${mut_path}
mkdir -p ${result_dir}/analysis

for coh in ${cohs};do
  
  #sparse and generate mutation informations from vcf files into integrated format
  pid_file=${pid_path}/${coh}.txt
  mut_dir=${mut_path}/${coh} 
  job_parser=`qsub -N ${coh}_parser -v coh=${coh},roi=${roi_file},pid_file=${pid_file},mut_dir=${mut_dir},analysis_scripts_dir=${analysis_scripts_dir},result_dir=${result_dir},output_path=${mut_path} -j oe -o ${output_dir} ${analysis_scripts_dir}/information_parser.sh -l walltime=1:00:00`

  #calculate score per functional element
  bed_file=${result_dir}/${coh}.bed
  snv_out=${result_dir}/analysis/${coh}_snv.bed
  cnv_out=${result_dir}/analysis/${coh}_cnv.bed
  sco_out=${result_dir}/analysis/${coh}_score.bed
  sco_out2=${result_dir}/analysis/${coh}_score2.bed
  out_put=${result_dir}/analysis/${coh}_result.bed
  #rm ${out_put}
  job_convert=`echo "python ${analysis_scripts_dir}/matrix_convert.py ${bed_file} ${snv_out} ${cnv_out}"|qsub -N ${coh} -j oe -o ${output_dir} -l walltime=00:15:00 -W depend="afterok:"${job_parser}`
  #job_convert=`echo "python ${analysis_scripts_dir}/matrix_convert.py ${bed_file} ${snv_out} ${cnv_out}"|qsub -N ${coh} -j oe -o ${output_dir} -l walltime=00:15:00 `
  #job_calc=`echo "Rscript ${analysis_scripts_dir}/scoring.r ${cnv_out} ${snv_out} ${sco_out}"|qsub -N ${coh}_scoring -j oe -o ${output_dir} -l walltime=5:00:00 -W depend="afterok:"${job_convert}`
  #job_calc=`echo "Rscript ${analysis_scripts_dir}/scoring_promoter.r ${snv_out} ${sco_out}"|qsub -N ${coh}_scoring -j oe -o ${output_dir} -l walltime=1:00:00 -W depend="afterok:"${job_convert}`
  
  #generate qq_plot graph and select best candidates from 10 turns of randomization
  ##job_qqplot=`echo "Rscript ${analysis_scripts_dir}/plot_qq.r ${coh} ${sco_out} ${result_dir}/analysis " |qsub -N ${coh}_qqplot -j oe -o ${output_dir} -l walltime=3:00:00 -W depend="afterok:"${job_calc}`
  
  #integrate candidates
  ##job_integration=`echo "cat ${result_dir}/analysis/qqplot*_${coh}.txt > ${result_dir}/analysis/candidates_${coh}.txt"|qsub -N ${coh}_integration -j oe -o ${output_dir} -l walltime=0:02:00 -W depend="afterok:"${job_qqplot}`
  ##job_ranking=`echo "python ${analysis_scripts_dir}/candidates_ranking.py ${result_dir} ${coh}"|qsub -N ${coh}_ranking -j oe -o ${output_dir} -l walltime=0:10:00 -W depend="afterok:"${job_integration}`
  #job_ranking=`echo "python ${analysis_scripts_dir}/candidates_ranking.py ${result_dir} ${coh}"|qsub -N ${coh}_ranking -j oe -o ${output_dir} -l walltime=0:10:00 `
  
  
  #sort -k 5 -gr ${sco_out} > ${sco_out2}
  #mv ${sco_out2} ${sco_out}

  
  done
