# fin = open("/abi/data/chen/promoter_project/im/f3.txt","r")
# cgcs = []
# for line in fin:
# 	line = line.strip("\n")
# 	cgcs.append(line.split()[0])
# print(cgcs)
# fin.close()

# info = ["none"]*len(cgcs)
# promoters = open("/home/hongc/PhD/promoters/code/promoter_hg.txt","r")
# for line in promoters:
# 	line = line.strip("\n")
# 	el = line.split("\t")
# 	pos = el[0]+"\t"+el[1]+"\t"+el[2]
# 	gene = el[3]
# 	for i in range(len(cgcs)):
# 		if cgcs[i]==gene:
# 			info[i] = pos
# promoters.close()
# print(info)

# info[0] = "chr\tstart\tend"

# fout1 = open("/abi/data/chen/promoter_project/im/p1f3.txt","w")
# for i in range(len(cgcs)):
# 	print >>fout1, info[i]+"\t"+cgcs[i]

# fout1.close()

import sys
import os
#from pybedtools import BedTool

# os.system("paste /abi/data/chen/promoter_project/im/p1f3.txt /abi/data/chen/promoter_project/im/f3.txt > /abi/data/chen/promoter_project/im/merge.txt")

# fin2 = open("/abi/data/chen/promoter_project/im/merge.txt","r")
# for line in fin2:
# 	line = line.strip("\n")
# 	el = line.split("\t")
# 	#pid = el[5]
# 	if el[5]!="pid":
# 		print el[5]
# 		os.system("cp /ibios/co02/chen/pcawg_2018/snv_mnv_del/"+ el[5] +".txt /abi/data/chen/promoter_project/im/vcfs")

# fin2.close()

# linux: cd /abi/data/chen/promoter_project/im/vcfs
# linux: for i in *; do awk '{print $0"\t"FILENAME}' $i > $i.bk;done
# linux: cat *.bk > snvs.bed
# linux: mv snvs.bed test/
# module load bedtools/2.16.2
# cut -f 1,2,3,4,5,6 /abi/data/chen/promoter_project/im/merge.txt > /abi/data/chen/promoter_project/im/info.txt
# add one # before first line of info.txt
# bedtools intersect -a /abi/data/chen/promoter_project/im/info.txt -b /abi/data/chen/promoter_project/im/vcfs/test/snvs.bed -wb > /abi/data/chen/promoter_project/im/merge_pos_f3.txt

fin3 = open("/abi/data/chen/promoter_project/im/merge_pos_f3.txt","r")
fout2 = open("/abi/data/chen/promoter_project/im/f3_result.txt","w")
for line in fin3:
	line = line.strip("\n")
	el = line.split("\t")
	if el[5] == el[12].split(".")[0]:
		print >>fout2,line
fin3.close()
fout2.close()

# cut -f 7,8,9,4,6,10,11 /abi/data/chen/promoter_project/im/f3_result.txt > f3_with_pos.txt