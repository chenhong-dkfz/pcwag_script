import sys

cgc = open("/ibios/co02/chen/casino_code/coh_pid/others/cgc_gene_jul_2016.csv","r")
cgc_gene = []
for line in cgc:
	line = line.strip("\n")
	if line not in cgc_gene:
		cgc_gene.append(line)

#ftest = open(sys.argv[1],"r")
ftest = open("/ibios/co02/chen/new_casino/reports/list200/"+sys.argv[1]+"_"+sys.argv[2]+".txt","r")

overlap = 0
lq = 0
ids = []

for line in ftest:
	line = line.strip("\n")
	els = line.split(";")
	lq += 1
	flag = 0
	for el in els:
		if el in cgc_gene:
			flag = 1
	if flag == 1:
		ids.append(lq)
	overlap += flag

print "******"
print "result for "+sys.argv[1]+" of "+sys.argv[2]+" is "+str(overlap)+"/"+str(lq)+"!"
#print "the ranks are "+" ".join(ids)
print ids
print "*****"