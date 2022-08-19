#pid_path=/home/hongc/CaSINo/coh_pid 
result_dir=/ibios/co02/chen/new_casino3/results/promoter_randoms
#coh='PACA-CA'
rg_path=/ibios/co02/chen/casino_code/coh_pid/regions_of_interests/tools
roi=${rg_path}/promoter.txt
roi_snv=${rg_path}/promoter_forsnv.txt
roi_cnv=${rg_path}/promoter_forcnv.txt
pid_path=/home/hongc/CaSINo/coh_pid
pid_file=${pid_path}/${coh}.txt
mut_path=/ibios/co02/chen/random
mut_dir=${mut_path}/${coh}
code_dir=/home/hongc/CaSINo/code
analysis_scripts_dir=${code_dir}/analysis_scripts

while read line; do

mkdir -p ${mut_dir}/${line}
mkdir -p ${result_dir}
output_path=${mut_path}

#only for snv and cnv files

#bedtools intersect -a ${roi} -b ${mut_dir}/${line}/snv.bed -wao  > ${mut_dir}/${line}/snv_tmp.bed
bedtools intersect -a ${roi} -b ${mut_dir}/${line}/snv-${turn}.bed -wao  > ${mut_dir}/${line}/snv_tmp${turn}.bed
#python ${analysis_scripts_dir}/merge_same_positions.py ${mut_dir}/${line}/snv_tmp.bed ${mut_dir}/${line}/snv2.bed 1
python ${analysis_scripts_dir}/merge_same_positions.py ${mut_dir}/${line}/snv_tmp${turn}.bed ${mut_dir}/${line}/snv2${turn}.bed 1
#cut -f 5 ${mut_dir}/${line}/snv2.bed >  ${mut_dir}/${line}/snv3.bed
#rm ${mut_dir}/${line}/snv2.bed  ${mut_dir}/${line}/snv_tmp.bed
#mv ${mut_dir}/${line}/snv3.bed ${mut_dir}/${line}/new_snv.bed

cut -f 5 ${mut_dir}/${line}/snv2${turn}.bed >  ${mut_dir}/${line}/snv3${turn}.bed
rm ${mut_dir}/${line}/snv2${turn}.bed  ${mut_dir}/${line}/snv_tmp${turn}.bed
mv ${mut_dir}/${line}/snv3${turn}.bed ${mut_dir}/${line}/new_snv${turn}.bed

#sed -i '1s/^/'${line}'\n/' ${mut_dir}/${line}/new_snv.bed 
sed -i '1s/^/'${line}'\n/' ${mut_dir}/${line}/new_snv${turn}.bed 

done<${pid_file}

#paste ${roi_snv} `ls  ${mut_dir}/*/new_snv.bed` > ${result_dir}/${coh}.bed

paste ${roi_snv} `ls  ${mut_dir}/*/new_snv${turn}.bed` > ${result_dir}/${coh}_${turn}.bed





