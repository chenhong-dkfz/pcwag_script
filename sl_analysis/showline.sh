snv_path=/ibios/co02/chen/pcawg_2018/snv_mnv_del
cnv_path=/ibios/co02/chen/pcawg_2018/cnv
sv_path=/ibios/co02/chen/pcawg_2018/sv
plist=/home/hongc/project_august/rep.txt

while read pid; do

wc -l ${snv_path}/${pid}.txt
wc -l ${cnv_path}/${pid}.txt
wc -l ${sv_path}/${pid}.txt

done <${plist}

