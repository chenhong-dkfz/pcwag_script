#pid_path=/home/hongc/CaSINo/coh_pid 
#result_dir=/ibios/co02/chen/new_casino/results
#coh='MB2'
#roi=${pid_path}/coding.bed
##roi_sorted=${pid_path}/gencode_s.bed
#mut_path=/ibios/co02/chen/new_casino/intermediate
#pid_path=/home/hongc/CaSINo/coh_pid
#pid_file=${pid_path}/${coh}.txt
#mut_dir=${mut_path}/${coh}
#code_dir=/home/hongc/CaSINo/code
#analysis_scripts_dir=${code_dir}/analysis_scripts

#sortBed -i ${roi} > ${roi_sorted}
#uniq ${roi_sorted} > ${roi_sorted}2

while read line; do

mkdir -p ${mut_dir}/${line}
mkdir -p ${result_dir}
#output_path=${mut_path}
#echo "works"
python ${analysis_scripts_dir}/parser_march.py ${line} ${coh} ${output_path}

#only for snv files

bedtools intersect -a ${roi} -b ${mut_dir}/${line}/snv.bed -wao  > ${mut_dir}/${line}/snv_tmp.bed
#head -n 1000 snv_tmp.bed > snv_tmp2.bed
python ${analysis_scripts_dir}/merge_same_positions.py ${mut_dir}/${line}/snv_tmp.bed ${mut_dir}/${line}/snv2.bed 1
#mv ${mut_dir}/${line}/snv_tmp.bed ${mut_dir}/${line}/snv_tmp_rsv.bed
cut -f 5 ${mut_dir}/${line}/snv2.bed >  ${mut_dir}/${line}/snv3.bed
#mv  ${mut_dir}/${line}/snv_tmp.bed  ${mut_dir}/${line}/snv.bed 
sed -i '1s/^/'${line}'\n/' ${mut_dir}/${line}/snv3.bed 
wc -l < ${mut_dir}/${line}/snv.bed > ${mut_dir}/${line}/count.txt
cat ${mut_dir}/${line}/snv3.bed ${mut_dir}/${line}/count.txt > ${mut_dir}/${line}/snv4.bed
rm ${mut_dir}/${line}/count.txt

done<${pid_file}


roi2=${result_dir}/temp_roi_${coh}.txt
bot=${result_dir}/temp_bot_${coh}.txt
echo -e '#chromo\tall\tall\tgene' > ${bot}
cat ${roi} ${bot} > ${roi2}
rm ${bot}

paste ${roi2} `ls  ${mut_dir}/*/snv4.bed` > ${result_dir}/${coh}.bed
rm ${roi2}
#bedtools intersect -a ${roi} -b ${pid_path}/reptiming.weg -wao > ${result_dir}/replication_timing.bed
#bedtools intersect -a ${roi} -b ${pid_path}/hic.bed -wao > ${result_dir}/hic.bed
#cut -f 1,2,3,4,9 hic.bed > hic_2.bed
#cut -f 1,2,3,4,8 replication_timing.bed  > replication_timing_2.bed
#python ${analysis_scripts_dir}/merge_same_positions.py ${result_dir}/replication_timing_2.bed ${result_dir}/replication_timing.bed 2
#python ${analysis_scripts_dir}/merge_same_positions.py ${result_dir}/hic_2.bed ${result_dir}/hic.bed 2
