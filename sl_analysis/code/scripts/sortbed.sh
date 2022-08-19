
pid_file="/ibios/co02/chen/pcawg_2018/pid_list/pids.tsv"
module load bedtools/2.23.0
source_dir="/ibios/co02/chen/pcawg_2018"
while read line; do
    #echo ${line} 
    sortBed -i ${source_dir}/snv_us/${line}.txt > ${source_dir}/snv_exon/${line}.txt
   
done<${pid_file}
    