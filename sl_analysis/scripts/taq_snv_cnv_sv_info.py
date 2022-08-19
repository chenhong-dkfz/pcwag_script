import os
import subprocess

fin = open("/home/hongc/project_august/sample_info_tumor_aliquot.txt")
fout = open("/home/hongc/project_august/sample_stats.txt","w")

for line in fin:
	if line[0] != "#":
		line = line.strip("\n")
		el = line.split("\t")
		ta = el[10]

		if os.path.isfile("/ibios/co02/chen/pcawg_2018/snv_mnv_del/"+ta+".txt"):
			output = subprocess.check_output("wc -l /ibios/co02/chen/pcawg_2018/snv_mnv_del/"+ta+".txt", shell=True)
			n1 = output.split()[0] 
		else:
			n1 = "0"
		if os.path.isfile("/ibios/co02/chen/pcawg_2018/cnv/"+ta+".txt"):
			output = subprocess.check_output("wc -l /ibios/co02/chen/pcawg_2018/cnv/"+ta+".txt", shell=True)
			n2 = output.split()[0]
		else:
			n2 = "0"
		if os.path.isfile("/ibios/co02/chen/pcawg_2018/sv/"+ta+".txt"):
			output = subprocess.check_output("wc -l /ibios/co02/chen/pcawg_2018/sv/"+ta+".txt", shell=True)
			n3 = output.split()[0]
		else:
			n3 = "0"
		print >> fout, ta + "\t" + n1 + "\t" +n2 + "\t" + n3

fin.close()
fout.close()
