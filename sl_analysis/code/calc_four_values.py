	# mutation levels: real biallelic level ( homozygous deletion and heterozygous deletion + matched snv) code: 5
	# likely biallelic level (2 snvs) code: 4
	# monollelic level (1 snv and 1 heterozygous deletion) code: 3
	# strong amplification (cn > 8) code: 2
	# normal amplification (cn > 4 and cn < 8) code: 1
	# non mutation code: 0
import os



def calc_four(a,b,typa,typb):
	if typa == "ra":
		a=a.replace("1","0")
		a=a.replace("2","0")
		a=a.replace("3","0")
		a=a.replace("4","0")
		a=a.replace("5","1")
	elif typa == "pa":
		a=a.replace("1","0")
		a=a.replace("2","0")
		a=a.replace("3","0")
		a=a.replace("4","1")
		a=a.replace("5","1")
	elif typa == "mut":
		a=a.replace("1","0")
		a=a.replace("2","0")
		a=a.replace("3","1")
		a=a.replace("4","1")
		a=a.replace("5","1")
	elif typa == "nam":
		a=a.replace("1","1")
		a=a.replace("2","1")
		a=a.replace("3","0")
		a=a.replace("4","0")
		a=a.replace("5","0")
	elif typa == "pam":
		a=a.replace("1","0")
		a=a.replace("2","1")
		a=a.replace("3","0")
		a=a.replace("4","0")
		a=a.replace("5","0")
	if typb == "ra":
		b=b.replace("1","0")
		b=b.replace("2","0")
		b=b.replace("3","0")
		b=b.replace("4","0")
		b=b.replace("5","1")
	elif typb == "pa":
		b=b.replace("1","0")
		b=b.replace("2","0")
		b=b.replace("3","0")
		b=b.replace("4","1")
		b=b.replace("5","1")
	elif typb == "mut":
		b=b.replace("1","0")
		b=b.replace("2","0")
		b=b.replace("3","1")
		b=b.replace("4","1")
		b=b.replace("5","1")
	elif typb == "nam":
		b=b.replace("1","1")
		b=b.replace("2","1")
		b=b.replace("3","0")
		b=b.replace("4","0")
		b=b.replace("5","0")
	elif typb == "pam":
		b=b.replace("1","0")
		b=b.replace("2","1")
		b=b.replace("3","0")
		b=b.replace("4","0")
		b=b.replace("5","0")

	a = a.split("\t")
	b = b.split("\t")
	full_length = len(a)
	hit_a = a.count("1")
	hit_b = b.count("1")
	hit_ab = 0
	for i in range(len(a)):
		if a[i] == "1" and b[i] == "1":
			hit_ab += 1
	return [str(full_length),str(hit_a),str(hit_b),str(hit_ab)]

#print(calc_four("2	2	2	3	1	0	0	1","0	1	2	1	1	1	0	0","pam","pam"))

cohorts = ["BLCA-US","CLLE-ES","GBM-US","LGG-US","MALY-DE","PAEN-AU","READ-US",
"BOCA-UK","CMDI-UK","HNSC-US","LICA-FR","MELA-AU","PBCA-DE","RECA-EU",
"BRCA-EU","COAD-US","KICH-US","LIHC-US","ORCA-IN","SARC-US",
"BRCA-UK","DLBC-US","KIRC-US","LINC-JP","OV-AU","SKCM-US",
"BRCA-US","EOPC-DE","KIRP-US","LIRI-JP","OV-US","PRAD-CA","STAD-US",
"BTCA-SG","ESAD-UK","LAML-KR","LUAD-US","PACA-AU","PRAD-UK","THCA-US",
"CESC-US","GACA-CN","LAML-US","LUSC-US","PACA-CA","PRAD-US","UCEC-US"]




datas = []
genes = [] # PTEN 24605th gene, CHD1L gene 2405th gene.

pid_coh = {}

ta_pro = open("/abi/data/chen/project_august/pcawg_ta_project.txt")
for line in ta_pro:
	if line[0:3]!= "tum":
		line = line.strip("\n")
		el = line.split("\t") 
		pid_coh[el[0]] = el[1]
ta_pro.close()



#fin = open("/abi/data/chen/project_august/mutation_pten.csv","r") # gene

fin = open("/abi/data/chen/project_august/mutation2.csv","r") # all
for line in fin:
	if line[0]=="g":
		line = line.strip("\n")
		pid_line = line.split("\t",1)[1]
		el = pid_line.split("\t")
		coh_line = []
		for el1 in el:
			coh_line.append(pid_coh[el1])


	if line[0]!="g":
		line = line.strip("\n")
		gene = line.split("\t",1)[0]
		genes.append(gene)
		line2 = line.split("\t",1)[1]
		datas.append(line2)

print genes[25077]  #TIAL1
#print genes[24604] PTEN
#print genes[2404]   # gene
#print genes[2404]

fin.close()

print len(genes)
print len(datas)

typs = ["nam","ra","pa","pam","mut"]

print len(datas)
print len(coh_line)

def data_redo(idx,cohort):
	string = datas[idx]
	string = string.split("\t")
	selected = [i for i, x in enumerate(coh_line) if x == cohort]
	string = [string[i] for i in selected]
	string = "\t".join(string)
	return(string)


#typ1 = "nam"
#typ2 = "ra"



for typ1 in typs:
	for typ2 in typs:
		print typ1+"_"+typ2
		for coh in cohorts:
			os.system("mkdir -p /abi/data/chen/project_august/mutation25_tial1/per_cohort/"+coh+"/TIAL1")
			fo = open("/abi/data/chen/project_august/mutation25_tial1/per_cohort/"+coh+"/TIAL1/"+typ1+"_"+typ2+".txt","w")   # gene
			for i in range(len(datas)):
				#if i == 3:
				#	print  datas[3]
				d1 = data_redo(25077,coh)   # gene
				#d1 = data_redo(2404,coh)   # gene
				#print(len(d1))
				#if i==3:
				#	print d1
				#	print "gogogo"
				d2 = data_redo(i,coh)
				ps = calc_four(d1,d2,typ1,typ2)     
				#ps = calc_four(datas[24604],datas[i],typ1,typ2)     # gene
				k = "\t".join(ps)
				print >> fo, k
			fo.close()

 
