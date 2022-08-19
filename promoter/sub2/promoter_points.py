# fcoding = open("/ibios/co02/chen/pcawg_2018/result/matrix/coding/select/mutation.csv","r")

# mut_gene = []

# for line in fcoding:
# 	el = line.split("\t")
# 	if el[0] != "gene_name":
# 		mut_gene.append(el[0])

# fcoding.close()

# exp_gene = []

# fexp = open("/abi/data/chen/promoter_project/exp_file_ta.csv","r")

# for line in fexp:
# 	el = line.split("\t")
# 	if el[0] != "#gene":
# 		exp_gene.append(el[0])
# fexp.close()

# i = 0
# j = 0

# for gene in mut_gene:
# 	if gene not in exp_gene:
# 		print gene
# 		i = i + 1


# for gene in exp_gene:
# 	if gene not in mut_gene:
# 		#print gene
# 		j = j + 1
# 		if "HA" not in gene:
# 			print gene

# print i
# print j

# frep = open("/abi/data/chen/promoter_project/replicated_genes.txt","w")

# for g in exp_gene:
# 	if g[-2:]=="HA" and g not in mut_gene:
# 		print >>frep,g[:-2]

# frep.close()


##################################################################################

f_coding = "/abi/data/chen/promoter_project/coding_amp.csv"
f_exp = "/abi/data/chen/promoter_project/exp_file_ta.csv"
f_pm = "/abi/data/chen/promoter_project/pm_snv.csv"
f_rg = "/abi/data/chen/promoter_project/replicated_genes.txt"
#f_gl = "/abi/data/chen/promoter_project/gene_list.txt"
f_gl = "/home/hongc/PhD/promoters/code/promoter_hg.txt"

gene_list = []
gl = open(f_gl,"r")
for line in gl:
	if line[0] != "#":
		line = line.strip("\n")
		gene0 = line.split("\t")[3]
		gene_list.append(gene0)
gl.close()

print gene_list[0]
print gene_list[1]

pm = open(f_pm,"r")
i = 0
pm_pid_list = []
pm_gene_list = []
for line in pm:
	if i==0:
		pid = []
		line = line.strip("\n")
		els = line.split("\t")
		for el in els:
			pid.append(el)
		i = i + 1
	else:
		line = line.strip("\n")
		els = line.split("\t")
		for j in range(len(els)):
			if els[j]=="1":
				pm_pid_list.append(pid[j])
				pm_gene_list.append(gene_list[i-1])
		i = i+1

pm.close()
#print pm_pid_list
print len(pm_pid_list)
print i

exp = open(f_exp,"r")
i = 0
exp_pid_list = []
exp_gene_list = []
exp_zscore_list = []
for line in exp:
	if i==0:
		pid = []
		line = line.strip("\n")
		els = line.split("\t")
		for el in els:
			pid.append(el)
		i = i + 1
	else:
		line = line.strip("\n")
		els = line.split("\t")
		for j in range(len(els)):
			if j != 0:
				if float(els[j]) > 3:
					exp_pid_list.append(pid[j])
					exp_gene_list.append(els[0])
					exp_zscore_list.append(els[j])
			else:
				pass
		i = i+1

exp.close()
#print exp_pid_list
print len(exp_pid_list)
print i


coding = open(f_coding,"r")
i = 0
coding_pid_list = []
coding_gene_list = []
#exp_zscore_list = []
for line in coding:
	if i==0:
		pid = []
		line = line.strip("\n")
		els = line.split("\t")
		for el in els:
			pid.append(el)
		i = i + 1
	else:
		line = line.strip("\n")
		els = line.split("\t")
		for j in range(len(els)):
			if els[j]=="1":
				coding_pid_list.append(pid[j])
				coding_gene_list.append(gene_list[i-1])
		i = i+1

coding.close()
#print exp_pid_list


print len(pm_pid_list)
print len(pm_gene_list)
print len(exp_pid_list)
print len(exp_zscore_list)
print len(exp_gene_list)
print len(coding_pid_list)
print len(coding_gene_list)

pm_combi = []
for i in range(len(pm_pid_list)):
	pm_combi.append(pm_pid_list[i]+"_"+pm_gene_list[i])

exp_combi = []
for i in range(len(exp_pid_list)):
	exp_combi.append(exp_pid_list[i]+"_"+exp_gene_list[i])

coding_combi = []
for i in range(len(coding_pid_list)):
	coding_combi.append(coding_pid_list[i]+"_"+coding_gene_list[i])

fout = open("/abi/data/chen/promoter_project/fout.txt","w")



intersect_combi = list(set(pm_combi).intersection(exp_combi))
print(len(intersect_combi))

exp_dict = {}

print "dict adding..."

for i in range(len(exp_combi)):
	exp_dict[exp_combi[i]] = exp_zscore_list[i]

print "dict printing..."

kk = 0
for j in intersect_combi:
	print kk
	kk = kk + 1
	if j in exp_dict.keys():
		print >>fout, j+","+exp_dict[j]


fout.close()


print "combi2 printing..."

intersect_combi2 = list(set(coding_combi).intersection(intersect_combi))
fout2 = open("/abi/data/chen/promoter_project/fout2.txt","w")
for i in intersect_combi2:
	print >>fout2,i

fout2.close()