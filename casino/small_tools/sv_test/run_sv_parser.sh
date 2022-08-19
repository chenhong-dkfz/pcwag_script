pid_path=/home/hongc/CaSINo/coh_pid
cohs="PBCA-DE PACA-AU PACA-CA PAEN-AU SARC-US OV-AU OV-US BRCA-EU BRCA-UK PRAD-UK PRAD-US EOPC-DE" 
#cohs="EOPC-DE"
#pid_file=${pid_path}/${coh}.txt
script_path=/home/hongc/CaSINo/code/small_tools/sv_test/sv_analysis.py
#result_path=/home/hongc/CaSINo/code/small_tools/sv_test/result/promoter
#coh='sv_pid_list'


roi=${pid_path}/enhancer.txt



for coh in ${cohs};do

pid_file=/home/hongc/CaSINo/coh_pid2/${coh}.txt
code_dir=/home/hongc/CaSINo/code
analysis_scripts_dir=${code_dir}/analysis_scripts

result_path=/home/hongc/CaSINo/code/small_tools/sv_test/result/enhancer/${coh}

mkdir -p ${result_path}



while read line;do
python ${script_path} ${line} ${result_path}/${line}.raw.bed

sortBed -i ${result_path}/${line}.raw.bed|uniq|mergeBed -nms -i > ${result_path}/${line}.sv_merged.bed
bedtools intersect -a ${roi} -b ${result_path}/${line}.sv_merged.bed -wao  > ${result_path}/${line}.mgd_tmp.bed
python ${analysis_scripts_dir}/merge_same_positions.py ${result_path}/${line}.mgd_tmp.bed ${result_path}/${line}.sv.bed 1
#rm ${result_path}/${line}.raw.bed
#rm ${result_path}/${line}.merged.bed
rm ${result_path}/${line}.sv_merged.bed
rm ${result_path}/${line}.mgd_tmp.bed

cut -f 5 ${result_path}/${line}.sv.bed > ${result_path}/${line}.sv2.bed
rm ${result_path}/${line}.sv.bed
find ${result_path} -size 0 -delete

if [ -f ${result_path}/${line}.sv2.bed ];then
sed -i '1s/^/'${line}'\n/' ${result_path}/${line}.sv2.bed
python /home/hongc/CaSINo/code/analysis_scripts/cnvsum.py ${result_path}/${line}.raw.bed ${result_path}/${line}.sumcnv.bed
cat ${result_path}/${line}.sv2.bed ${result_path}/${line}.sumcnv.bed > ${result_path}/${line}.sv3.bed
fi

rm ${result_path}/${line}.raw.bed
rm ${result_path}/${line}.sv2.bed
rm ${result_path}/${line}.sumcnv.bed

done<${pid_file}

find ${result_path} -size 0 -delete
find ${result_path} -size 2 -delete
roi1=/home/hongc/CaSINo/code/small_tools/sv_test/temp/temp_roi1_${coh}.txt
cp ${roi} ${roi1}
roi2=/home/hongc/CaSINo/code/small_tools/sv_test/temp/temp_roi_${coh}.txt
bot=/home/hongc/CaSINo/code/small_tools/sv_test/temp/temp_bot_${coh}.txt
echo -e '#chromo\tall\tall\tgene' > ${bot}
#sed -i '1s/^/'#chromo\tall\tall\tgene'\n/' ${roi1}
cat ${bot} ${roi1} ${bot} > ${roi2}
rm ${bot}
rm ${roi1}

paste ${roi2} `ls  ${result_path}/*.sv3.bed` > ${result_path}/${coh}.bed
rm ${roi2}
rm ${result_path}/*.sv3.bed
done
