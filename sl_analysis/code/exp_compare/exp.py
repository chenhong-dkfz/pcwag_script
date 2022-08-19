gene1 = "PTEN"
gene2 = "CHD1"

gc19 = open("/ibios/co02/chen/pcawg_2018/regions/gencode.v19.txt","r")
for line in gc19:
	line = line.strip("\n")
	el = line.split("\t")
	if el[1] == gene1:
		ensg1 = el[0]
	if el[1] == gene2:
		ensg2 = el[0]

gc19.close()

gene_list = [gene1,gene2]

zsc = open("/abi/data/chen/expression_final_version/fpkm/fpkm.csv","r")
fo = open("/abi/data/chen/project_august/exp_comp/im.txt","w")

for line in zsc:
	if line[0] == "#":
		line = line.strip("\n")
		print >>fo,line
	else:
		line = line.strip("\n")
		line_gene = line.split("\t")[0]
		if line_gene == ensg1:
			#print >>fk,line
			line = line.replace(line_gene,gene1)
			print >>fo,line
		if line_gene == ensg2:
			line = line.replace(line_gene,gene2)
			print >>fo,line

zsc.close()
fo.close()


