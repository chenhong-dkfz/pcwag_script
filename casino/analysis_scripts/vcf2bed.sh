#pid_path=/home/hongc/CaSINo/coh_pid 
result_dir=/ibios/co02/chen/new_casino/results/promoter
coh='BLCA-US SKCM-US THCA-US'
rg_path=/ibios/co02/chen/casino_code/coh_pid/regions_of_interests/tools
roi=${rg_path}/promoter.txt
roi_snv=${rg_path}/promoter_forsnv.txt
roi_cnv=${rg_path}/promoter_forcnv.txt
pid_path=/home/hongc/CaSINo/coh_pid
pid_file=${pid_path}/${coh}.txt
mut_path=/ibios/co02/chen/data2015_3new
mut_dir=${mut_path}/${coh}
code_dir=/home/hongc/CaSINo/code
analysis_scripts_dir=${code_dir}/analysis_scripts

while read line; do

mkdir -p ${mut_dir}/${line}
mkdir -p ${result_dir}
output_path=${mut_path}
#echo "works"
python ${analysis_scripts_dir}/parser.py ${line} ${coh} ${output_path}

#only for snv and cnv files

bedtools intersect -a ${roi} -b ${mut_dir}/${line}/snv_mnv_test_pass.bed -wao  > ${mut_dir}/${line}/snv_tmp.bed
python ${analysis_scripts_dir}/merge_same_positions.py ${mut_dir}/${line}/snv_tmp.bed ${mut_dir}/${line}/snv2.bed 1
cut -f 5 ${mut_dir}/${line}/snv2.bed >  ${mut_dir}/${line}/snv3.bed
#sed -i '1s/^/'${line}'\n/' ${mut_dir}/${line}/snv3.bed 
#wc -l < ${mut_dir}/${line}/snv_mnv_test_pass.bed > ${mut_dir}/${line}/count.txt
#cat ${mut_dir}/${line}/snv3.bed ${mut_dir}/${line}/count.txt > ${mut_dir}/${line}/snv4.bed
rm ${mut_dir}/${line}/snv_mnv_test_pass.bed ${mut_dir}/${line}/snv2.bed 
mv ${mut_dir}/${line}/snv3.bed ${mut_dir}/${line}/snv.bed

bedtools intersect -a ${roi} -b ${mut_dir}/${line}/cnv_test_pass.bed -wao > ${mut_dir}/${line}/cnv_tmp.bed
python ${analysis_scripts_dir}/merge_same_positions.py ${mut_dir}/${line}/cnv_tmp.bed ${mut_dir}/${line}/cnv.bed 1
cut -f 5 ${mut_dir}/${line}/cnv.bed >  ${mut_dir}/${line}/cnv_tmp.bed
mv  ${mut_dir}/${line}/cnv_tmp.bed  ${mut_dir}/${line}/cnv.bed
python ${analysis_scripts_dir}/cnvsum.py ${mut_dir}/${line}/cnv_test_pass.bed ${mut_dir}/${line}/sumcnv.bed
rm ${mut_dir}/${line}/cnv_test_pass.bed

#paste -d ';' ${mut_dir}/${line}/snv.bed ${mut_dir}/${line}/indel.bed ${mut_dir}/${line}/sv.bed ${mut_dir}/${line}/cnv.bed > ${mut_dir}/${line}/roi.bed
paste -d ';' ${mut_dir}/${line}/snv.bed ${mut_dir}/${line}/cnv.bed > ${mut_dir}/${line}/roi.bed
sed -i -- 's/\(;;\|;;;\)/;/g' ${mut_dir}/${line}/roi.bed
sed -i -- 's/\(^;\|;$\)//g' ${mut_dir}/${line}/roi.bed
sed -i -- 's/\(\.;\|;\.\)//g' ${mut_dir}/${line}/roi.bed
sed -i '1s/^/'${line}'\n/' ${mut_dir}/${line}/roi.bed 
cat ${mut_dir}/${line}/roi.bed ${mut_dir}/${line}/sumcnv.bed > ${mut_dir}/${line}/roi2.bed 

done<${pid_file}

paste ${roi_cnv} `ls  ${mut_dir}/*/roi2.bed` > ${result_dir}/${coh}.bed






