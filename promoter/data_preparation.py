import os
import sys

# gene list
add_a_gene = sys.argv[1]

target_genes = ["CDK12"]  # hugo symbols
target_genes.append(add_a_gene)

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

# gene_id_convert

gene_id = {} # hugo as key and ensg as value
gene_id_r = {} # ensg as key

gid = open("/ibios/co02/chen/pcawg_2018/regions/gencode.v19.txt","r")
for line in gid:
	line = line.strip("\n")
	el = line.split("\t")
	if el[1] not in gene_id.keys():
		gene_id[el[1]] = el[0]
	else:
		pass
gid.close()	

ensg_list = []

for target in target_genes:
	if target in gene_id.keys():
		gene_ensg = gene_id[target]
		ensg_list.append(gene_ensg)
		gene_id_r[gene_ensg] =  target
print ensg_list[1]

zsc = open("/abi/data/chen/expression_final_version/zscore/zscore.csv","r")
fo = open("/ibios/co02/chen/pcawg_2018/result/file_cache/exp_file_star.csv","w")
for line in zsc:
	if line[0] == "#":
		line = line.strip("\n")
		print >>fo,line
	else:
		line = line.strip("\n")
		line_gene = line.split("\t")[0]
		if line_gene in ensg_list:
			line = line.replace(line_gene,gene_id_r[line_gene])
			print >>fo,line

zsc.close()
fo.close()

fin = open("/ibios/co02/chen/pcawg_2018/result/file_cache/exp_file_star.csv","r")
fout = open("/ibios/co02/chen/pcawg_2018/result/file_cache/exp_file_ta.csv","w")

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

cohort_list = ["BLCA-US","CESC-US","ESAD-UK","KIRP-US","LIRI-JP","OV-US","PRAD-US","UCEC-US","BOCA-UK","CLLE-ES",
				"GACA-CN","LAML-KR","LUAD-US","PACA-AU","READ-US","BRCA-EU","CMDI-UK","GBM-US","LAML-US","LUSC-US",
				"PACA-CA","SARC-US","BRCA-UK","COAD-US","HNSC-US","LGG-US","MALY-DE","PAEN-AU","SKCM-US","BRCA-US",
				"DLBC-US","KICH-US","LICA-FR","ORCA-IN","PBCA-DE","STAD-US","BTCA-SG","EOPC-DE","KIRC-US","LIHC-US","OV-AU",
				"PRAD-UK","THCA-US","LINC-JP","PRAD-CA","RECA-EU","MELA-AU"]

# grenerate mutation matrix

fin = open("/ibios/co02/chen/pcawg_2018/result/matrix/coding/select/mutation.csv","r")
fo = open("/ibios/co02/chen/pcawg_2018/result/file_cache/mut1.csv","w")
for line in fin:
	if line[0] == "g":
		line = line.strip("\n")
		line = "#"+line
		print >>fo,line
	else:
		line = 	line.strip("\n")
		el = line.split("\t")
		if el[0] in target_genes:
			print >>fo,line
fin.close()
fo.close()

fin = open("/ibios/co02/chen/pcawg_2018/result/file_cache/mut1.csv","r")
fout = open("/ibios/co02/chen/pcawg_2018/result/file_cache/mut.csv","w")

for line in fin:
	if line[0] == "#":
		line = line.strip("\n")
		print >>fout, line
	else:
		line = line.strip("\n")
		# comment for filter
		line = line.replace("CNV_GAIN_1_","CNV_SG_")
		line = line.replace("CNV_GAIN_2_","CNV_SG_")
		line = line.replace("CNV_GAIN_3_","CNV_SG_")
		line = line.replace("CNV_GAIN_4_","CNV_SG_")
		# comment below for cnv5-8
		line = line.replace("CNV_GAIN_5_","CNV_SG_")
		line = line.replace("CNV_GAIN_6_","CNV_SG_")
		line = line.replace("CNV_GAIN_7_","CNV_SG_")
		line = line.replace("CNV_GAIN_8_","CNV_SG_")
		print >>fout,line
fin.close()
fout.close()




os.remove("/ibios/co02/chen/pcawg_2018/result/file_cache/mut1.csv")