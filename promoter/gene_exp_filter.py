fin = open("/abi/data/chen/expression_final_version/fpkm/fpkm.csv","r")
fo = open("/ibios/co02/chen/pcawg_2018/result/matrix/promoter/report/median_exp.txt","w")

import numpy

for line in fin:
	if line[0] != "#":
		line = line.strip("\n")
		el = line.split("\t")
		exp = []
		for i in range(len(el)):
			if i != 0:
				exp.append(float(el[i]))
		exp_median = numpy.median(exp)
		print >>fo,el[0]+"\t"+str(exp_median)

fin.close()
fo.close()