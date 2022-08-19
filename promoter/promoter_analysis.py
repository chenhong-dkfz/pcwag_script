gene_mapping = {}
start0 = []
gene0 = []

fin1 = open("/ibios/co02/chen/pcawg_2018/result/matrix/promoter/report/stats.txt","r")
i = 0
for line in fin1:
	print line
	if i < 47979:
		if line[0] != "#":
			line = line.strip("\n")
			el = line.split("\t")
			gene_mapping[el[5]] = el[3]
			i += 1
			start0.append(el[1])
			gene0.append(el[3])
	else:
		if "TERT\t" in line or "MYB\t" in line:
			print line
			line = line.strip("\n")
			el = line.split("\t")
			gene_mapping[el[5]] = el[3]
			i += 1
			start0.append(el[1])
			gene0.append(el[3])		
	

fin1.close()
print(len(gene_mapping))



# pid_convert

match = {} # tumor aliquote as value and star id as key
f_pid_rna = open("/ibios/co02/chen/pcawg_2018/regions/tumor_rna_new_uniq.txt","r")
for line in f_pid_rna:
	if line[0] != "T":
		line = line.strip("\n")
		el = line.split("\t")
		if el[1] not in match.keys():
			if el[1] != "":
				match[el[1]] = el[0]
		else:
			pass	
f_pid_rna.close()



zsc = open("/abi/data/chen/expression_final_version/zscore/zscore.csv","r")
fo = open("/ibios/co02/chen/pcawg_2018/result/matrix/promoter/report/exp_file_star.csv","w")
fk = open("/ibios/co02/chen/pcawg_2018/result/matrix/promoter/report/exp_file_fk.csv","w")

for line in zsc:
	if line[0] == "#":
		line = line.strip("\n")
		print >>fo,line
	else:
		line = line.strip("\n")
		line_gene = line.split("\t")[0]
		if line_gene in gene_mapping.keys():
			print >>fk,line
			line = line.replace(line_gene,gene_mapping[line_gene])
			print >>fo,line

zsc.close()
fo.close()
fk.close()


fs = open("/ibios/co02/chen/pcawg_2018/result/matrix/promoter/report/exp_file_fk.csv","w")
for x in gene_mapping.keys():
	print >>fs,x
	print >>fs,gene_mapping[x]
fs.close()

for i in match.keys():
	if "," in match[i]:
		del match[i]


print(len(match))

pm = open("/ibios/co02/chen/pcawg_2018/result/matrix/promoter/report/full_pm.csv","r")
pmo = open("/ibios/co02/chen/pcawg_2018/result/matrix/promoter/report/pm_top.csv","w")


for line in pm:
	if line[0] == "#":
		line = line.strip("\n")
		print >>pmo,line
	else:
		line = line.strip("\n")
		el = line.split("\t")
		for j in gene0:
			if j == el[3]:
				k = gene0.index(j)
				if start0[k] == el[1]:
					print >>pmo,line

pm.close()
pmo.close() 










fin = open("/ibios/co02/chen/pcawg_2018/result/matrix/promoter/report/exp_file_star.csv","r")
fout = open("/ibios/co02/chen/pcawg_2018/result/matrix/promoter/report/exp_file_ta.csv","w")

for line in fin:
	if line[0] == "#":
		line = line.strip("\n")
		els = line.split("\t")
		elo = []
		for el in els:
			if el == "#gene":
				elo.append("#gene")
			elif el in match.keys():
				elo.append(match[el])
			else:
				elo.append("known")
		line2 = "\t".join(elo)
		print >>fout,line2
	else:
		line = line.strip("\n")
		print >>fout,line

fin.close()
fout.close()

# muts = open("/ibios/co02/chen/pcawg_2018/result/matrix/promoter/report/pm.csv")

# with open("/ibios/co02/chen/pcawg_2018/result/matrix/promoter/report/pm.csv") as f:
#     first_line = f.readline()
# f.close()


