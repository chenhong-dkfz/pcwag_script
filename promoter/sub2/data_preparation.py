fin =  open("/ibios/co02/chen/pcawg_2018/result/matrix/coding/select/mutations.csv","r")
fout =  open("/ibios/co02/chen/pcawg_2018/result/matrix/coding/select/coding_amp.csv","w")


for line in fin:
	if line[0] == "0":
		line = line.strip("\n")
		print >>fout,line
	else:
		line = line.strip("\n")
		el = line.split("\t")
		new_line = []
		for els in el:
			if "CNV_GAIN" in els:
				new_line.append("1")
			else:
				new_line.append("0")
		newline = "\t".join(new_line)
		print >>fout,newline

fin.close()
fout.close()