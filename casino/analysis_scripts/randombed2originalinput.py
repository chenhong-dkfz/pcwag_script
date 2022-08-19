import os
import sys

############
#parameters#
############
sel_coh = sys.argv[1]
turn = sys.argv[2]

coh_pid_map = open("/home/hongc/project/codes/random_qq/coh_pid_mapping.csv","r")
#pid_list = open("/home/hongc/project/sample/coh_pid/"+coh+".txt","r")
pid_don = {}
for line in coh_pid_map:
	line = line.strip("\n")
	#print line
	elements = line.split("\t")
	coh = elements[0]
	doner_id = elements[1]
	pid = elements[2]
	pid_don[doner_id] = pid
"""
	try:
		os.mkdir("/ibios/co02/chen/random/"+coh)
	except:
		pass
	try:
		os.mkdir("/ibios/co02/chen/random/"+coh+"/"+pid)
	except:
		pass
"""

print "1st step done!"

#mut_path = "/icgc/dkfzlsdf/analysis/B080/crg/PanCan/data/output/annotated_set/SNV_Rand_U-D5000/"+sel_coh+"/"+sel_coh+".t-1.bed"
mut_path = "/icgc/dkfzlsdf/analysis/B080/crg/PanCan/data/output/2015-10-08-SantaCruz_withSupplement/complete_set/SNV_Rand_U-D25000/"+sel_coh+"/"+sel_coh+".t-"+turn+".bed"
#output_path = "/ibios/co02/chen/random/"+sel_coh+"/"+pid+"/"+"snv-"+turn+".bed"
#if os.path.isfile(output_path):
#	os.remove(output_path)

patient_info = {}

print "old file removed!"

mut = open(mut_path,"r")
for line in mut:
	if line[0]!= "C":
		line = line.strip("\n")
		elements = line.split("\t")
		chrom = elements[0]
		start = elements[1]
		end = elements[2]
		doner_id = elements[8]
		inf = chrom+"\t"+start+"\t"+str(int(end)+1)+"\tSNV"

		if doner_id in pid_don.keys():
			#print "doing!"
			pid = pid_don[doner_id]
			#print pid
			os.system("mkdir -p /ibios/co02/chen/random/"+sel_coh+"/"+pid)
			
			if pid not in patient_info.keys():
				patient_info[pid]=[]
				patient_info[pid].append(inf)
			else:
				patient_info[pid].append(inf)

			#print output
			#print >>output,chrom+"\t"+start+"\t"+str(int(end)+1)+"\tSNV"
			#output.close()
		else:
			pass
			#print doner_id

print "2nd step done!"

for patient in patient_info.keys():
	os.system("mkdir -p "+"/ibios/co02/chen/random/"+sel_coh+"/"+patient)
	output_path = open("/ibios/co02/chen/random/"+sel_coh+"/"+patient+"/"+"snv-"+turn+".bed","w")
	#output = open("/ibios/co02/chen/random/"+sel_coh+"/"+pid+"/"+"snv.bed","a")
	for l in patient_info[patient]:
		print >>output_path,l
	output_path.close()
	
print "written done!"