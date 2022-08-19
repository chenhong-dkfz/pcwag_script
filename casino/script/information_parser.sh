#parameters given
#roi=/ibios/co02/chen/tmp/roadmap_stringent_rois.bed.txt #positions of rois
#coh='LICA-FR' #cohort
#pid_file="/home/hongc/project/sample/coh_pid/${coh}.txt" #pids per cohort
#mut_dir=/ibios/co02/chen/roi/${coh}
#result_dir=/ibios/co02/chen/roi/results
#line='553912c0-8edf-4db2-a622-95a081a53519'
#analysis_scripts_dir=${code_dir}/analysis_scripts

while read line; do

mkdir -p ${mut_dir}/${line}
mkdir -p ${result_dir}
#parser informations from vcf files
python ${analysis_scripts_dir}/parser.py ${line} ${coh} ${output_path}

#merge roi with snv file
bedtools intersect -a ${roi} -b ${mut_dir}/${line}/snv_mnv_test_pass.bed -wao > ${mut_dir}/${line}/snv_tmp.bed
python ${analysis_scripts_dir}/merge_same_positions2.py ${mut_dir}/${line}/snv_tmp.bed ${mut_dir}/${line}/snv.bed
mv ${mut_dir}/${line}/snv_tmp.bed ${mut_dir}/${line}/snv_tmp_rsv.bed
cut -f 5 ${mut_dir}/${line}/snv.bed >  ${mut_dir}/${line}/snv_tmp.bed
mv  ${mut_dir}/${line}/snv_tmp.bed  ${mut_dir}/${line}/snv.bed

##merge roi with indel file 
bedtools intersect -a ${roi} -b ${mut_dir}/${line}/indel_test_pass.bed -wao > ${mut_dir}/${line}/indel_tmp.bed

cut -f 1,2,3,8 ${mut_dir}/${line}/indel_tmp.bed > ${mut_dir}/${line}/indel_tmp2.bed
uniq ${mut_dir}/${line}/indel_tmp2.bed > ${mut_dir}/${line}/indel.bed
mv ${mut_dir}/${line}/indel_tmp.bed ${mut_dir}/${line}/indel_tmp_rsv.bed
cut -f 4 ${mut_dir}/${line}/indel.bed >  ${mut_dir}/${line}/indel_tmp.bed
mv  ${mut_dir}/${line}/indel_tmp.bed  ${mut_dir}/${line}/indel.bed
sed -i -- 's/\.//g' ${mut_dir}/${line}/indel.bed
rm ${mut_dir}/${line}/indel_tmp2.bed

##merge roi with sv file
bedtools intersect -a ${roi} -b ${mut_dir}/${line}/sv_test_pass.bed -wao > ${mut_dir}/${line}/sv_tmp.bed

cut -f 1,2,3,8 ${mut_dir}/${line}/sv_tmp.bed > ${mut_dir}/${line}/sv_tmp2.bed
uniq ${mut_dir}/${line}/sv_tmp2.bed > ${mut_dir}/${line}/sv.bed
mv ${mut_dir}/${line}/sv_tmp.bed ${mut_dir}/${line}/sv_tmp_rsv.bed
cut -f 4 ${mut_dir}/${line}/sv.bed >  ${mut_dir}/${line}/sv_tmp.bed
mv  ${mut_dir}/${line}/sv_tmp.bed  ${mut_dir}/${line}/sv.bed
sed -i -- 's/\.//g' ${mut_dir}/${line}/sv.bed
rm ${mut_dir}/${line}/sv_tmp2.bed

##merge roi with cnv file
bedtools intersect -a ${roi} -b ${mut_dir}/${line}/cnv_test_pass.bed -wao > ${mut_dir}/${line}/cnv_tmp.bed
cp ${mut_dir}/${line}/cnv_tmp.bed ${mut_dir}/${line}/cnv_11.bed
python ${analysis_scripts_dir}/merge_same_positions2.py ${mut_dir}/${line}/cnv_tmp.bed ${mut_dir}/${line}/cnv.bed
cut -f 5 ${mut_dir}/${line}/cnv.bed >  ${mut_dir}/${line}/cnv_tmp.bed
#mv ${mut_dir}/${line}/cnv.bed ${mut_dir}/${line}/cnv_rsv.bed
mv  ${mut_dir}/${line}/cnv_tmp.bed  ${mut_dir}/${line}/cnv.bed
python ${analysis_scripts_dir}/cnvsum.py ${mut_dir}/${line}/cnv_test_pass.bed ${mut_dir}/${line}/sumcnv.bed

##paste all four files together and add a hat
paste -d ';' ${mut_dir}/${line}/snv.bed ${mut_dir}/${line}/indel.bed ${mut_dir}/${line}/sv.bed ${mut_dir}/${line}/cnv.bed > ${mut_dir}/${line}/roi.bed
sed -i -- 's/\(;;\|;;;\)/;/g' ${mut_dir}/${line}/roi.bed
sed -i -- 's/\(^;\|;$\)//g' ${mut_dir}/${line}/roi.bed
sed -i -- 's/\(\.;\|;\.\)//g' ${mut_dir}/${line}/roi.bed
sed -i '1s/^/'${line}'\n/' ${mut_dir}/${line}/roi.bed 
cat ${mut_dir}/${line}/roi.bed ${mut_dir}/${line}/sumcnv.bed > ${mut_dir}/${line}/roi2.bed 
done<${pid_file}


paste ${roi} `ls  ${mut_dir}/*/roi2.bed` > ${result_dir}/${coh}.bed
