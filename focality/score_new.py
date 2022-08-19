cohort_list = ["BLCA-US","CESC-US","ESAD-UK","KIRP-US","LIRI-JP","OV-US","PRAD-US","UCEC-US","BOCA-UK","CLLE-ES",
				"GACA-CN","LAML-KR","LUAD-US","PACA-AU","READ-US","BRCA-EU","CMDI-UK","GBM-US","LAML-US","LUSC-US",
				"PACA-CA","SARC-US","BRCA-UK","COAD-US","HNSC-US","LGG-US","MALY-DE","PAEN-AU","SKCM-US","BRCA-US",
				"DLBC-US","KICH-US","LICA-FR","ORCA-IN","PBCA-DE","STAD-US","BTCA-SG","EOPC-DE","KIRC-US","LIHC-US","OV-AU",
				"PRAD-UK","THCA-US","LINC-JP","PRAD-CA","RECA-EU","MELA-AU"]

#cohort_list = ["BLCA-US"]

matrix_path = "/abi/data/chen/co02/pcawg_2018/result/matrix/coding/"
write_path = "/abi/data/chen/newhome/promoter_focallity/data/"

import numpy as np

def focal_score_1(nlist):
	if sum(nlist) == 0:
		fs = 0
	else:
		fs = 0
		for x in nlist:
			if x!= 0:
				fs_a = np.log10(max(nlist))-np.log10(x)
				fs += fs_a
	return fs



# for cohort in cohort_list:
# 	file_name = matrix_path+cohort+"2.bed"
# 	fin = open(file_name,"r")
# 	losses = []
# 	for line in fin:
# 		if line[0]!="#":
# 			line = line.strip("\n")
# 			els = line.split("\t")
# 			gene_loss = []
# 			for el in els:
# 				short_gain = 0
# 				short_loss = 0
# 				muts = el.split(";")
# 				for mut in muts:
# 					if "CNV_GAIN" in mut:
# 						pass
# 					if "CNV_LOSS" in mut:
# 						loss_type = mut.split("_")[2]
# 						loss_range = mut.split("_")[3]
# 						loss_start = loss_range.split(",")[0]
# 						loss_end = loss_range.split(",")[1]
# 						loss_length = float(loss_end)-float(loss_start)
# 						if short_loss == 0:
# 							short_loss = loss_length
# 						else:
# 							if loss_length < short_loss:
# 								short_loss = loss_length
# 				gene_loss.append(short_loss)					
# 			losses.append(gene_loss)
# 	fin.close()

# 	file_name2 = write_path+cohort+".txt"
# 	fout = open(file_name2,"w")
# 	print>>fout,cohort
# 	for x in losses:
# 		fs = focal_score_1(x)
# 		print >>fout,fs


# 	fout.close()

focal_threshold = 10000000
#focal_threshold = 10000000

file_name = "/abi/data/chen/co02/pcawg_2018/result/matrix/coding/select/mutation.csv"
fin = open(file_name,"r")
losses = []
counts = []
for line in fin:
	if line[0]!="#":
		line = line.strip("\n")
		els = line.split("\t")
		gene_loss = []
		for el in els:
			short_gain = 0
			short_loss = 0
			muts = el.split(";")
			for mut in muts:
				if "CNV_GAIN" in mut:
					pass
				if "CNV_LOSS" in mut:
					loss_type = mut.split("_")[2]
					loss_range = mut.split("_")[3]
					loss_start = loss_range.split(",")[0]
					loss_end = loss_range.split(",")[1]
					loss_length = float(loss_end)-float(loss_start)
					if short_loss == 0:
						short_loss = loss_length
					else:
						if loss_length < short_loss:
							short_loss = loss_length
			if short_loss <= focal_threshold:
				gene_loss.append(short_loss)					
		losses.append(gene_loss)
		gene_loss_count = sum(1 for i in gene_loss if i != 0)
		counts.append(gene_loss_count)
#print(len(losses))
#print(losses)
fin.close()

file_name2 = write_path+"all_fscore_1e7.txt"
fout = open(file_name2,"w")
print>>fout,"all"
for x in losses:
	fs = focal_score_1(x)
	print >>fout,fs
fout.close()

print("step 1 finished!")

file_name3 = write_path+"all_counts_1e7.txt"
fout = open(file_name3,"w")
print>>fout,"all"
for x in counts:
	fs = x
	print >>fout,fs
fout.close()
