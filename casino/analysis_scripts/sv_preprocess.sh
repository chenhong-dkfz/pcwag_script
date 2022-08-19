cohs="EOPC-DE"
roi=/home/hongc/CaSINo/coh_pid/promoter.txt
#dir=/home/hongc/CaSINo/code/small_tools/sv_test/result/promoter
output_file=/home/hongc/CaSINo/code/small_tools/sv_test/analysis/output.bed

if [ -f ${output_file} ]
 then
  rm ${output_file}
fi

dir2=/home/hongc/CaSINo/code/small_tools/sv_test/analysis
mkdir -p ${dir2}

for coh in ${cohs};do

dir=/home/hongc/CaSINo/code/small_tools/sv_test/result/promoter/${coh}
files=`ls ${dir}`
  for file in ${files};do
    cut -f 5 ${dir}/${file} > ${dir2}/${file}.tmp
  done


paste ${roi} `ls ${dir2}/*.tmp` > ${output_file} 

done

