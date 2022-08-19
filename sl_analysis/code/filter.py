muts = ["pam","ram","pa","ra","mut"]

match = []
for i in muts:
	for j in muts:
		match.append(i+"_"+j)

print(match)

fout = open("/abi/data/chen/project_august/result_1311/gene_list_alc1.txt","w")

for m in match:
	try:
		print >>fout, m
		fin = open("/abi/data/chen/project_august/mutation25/cohortwise/result/ALC1/"+m+".txt","r")
		for line in fin:
			if line[0] != "c":
				line = line.strip("\n")
				el = line.split("\t")
				flag = 0 
				if el[2] != "TRUE":
					flag = 1
				if float(el[3]) >= 0.05:
					flag = 1
				if float(el[4]) > 0.5:
					flag = 1
				if float(el[8]) > 3:
					flag = 1
				if flag == 0:
					print >>fout,line
		fin.close()
	except:
		pass
fout.close()