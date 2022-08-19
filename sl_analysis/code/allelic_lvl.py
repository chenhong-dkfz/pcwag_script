import re 


def dh_detect(strings): #return a list of hom,dh,sh,null
	strings = strings.replace("SNV_0","NON") # remove intragenic,intronic...
	strings = strings.replace("SNV_1","NON") # remove RNA intronic
	strings = strings.replace("SNV_2","NON") # remove ncRNA exonic
	strings = strings.replace("SNV_3","NON") # remove 3' and 5' utr snvs
	strings = strings.replace("INDEL_2","SNV_6") # non-frameshift indels
	strings = strings.replace("INDEL_3","SNV_6") # frameshift indels

	mutations = strings.split(";")

	# mutation levels: real biallelic level ( homozygous deletion and heterozygous deletion + matched snv) code: 5
	# likely biallelic level (2 snvs) code: 4
	# monollelic level (1 snv and 1 heterozygous deletion) code: 3
	# strong amplification (cn > 8) code: 2
	# normal amplification (cn > 4 and cn < 8) code: 1
	# non mutation code: 0


	if len(mutations)==0: # if no mutations at all, null + 1
		#return [0,0,0,1]
		return 0
	elif re.match(".*CNV_GAIN.*",strings):
		amps = re.findall("CNV_GAIN_\d+",strings)
		floorcn = 999
		for amp in amps:
			copy_number = int(re.match("CNV_GAIN_(\d+)",amp).groups()[0])
			if copy_number < floorcn:
				floorcn = copy_number
		if floorcn > 8:
			return 2
		elif floorcn > 4:
			return 1
		else:
			return 0

	elif re.match(".*CNV_LOSS_HOM.*",strings): # one cnv loss hom, hom,dh +1
		#return [1,1,0,0]
		return 5
	elif strings.count("CNV_LOSS_HET")+strings.count("SNV")==0: # no deletions and no snvs, null + 1
		#return [0,0,0,1]
		return 0
	elif strings.count("CNV_LOSS_HET") == 0 and strings.count("SNV") > 1: # no deletion but multiple snv
		#return [0,0,1,0]
		return 4
	elif strings.count("CNV_LOSS_HET") == 0 and strings.count("SNV") == 1: # no deletion but single snv
		#return [0,0,1,0]
		return 3
	elif strings.count("CNV_LOSS_HET") != 0 and strings.count("SNV") == 0: # no snv, sh + 1		
		#return [0,0,1,0]
		return 3
	else: # at least one snv and one deletion
		#return 5

		deletions = []
		SNV = []
		for m in mutations:
			if m.count("CNV_LOSS_HET") > 0:
				deletions.append(m.split("_")[3])
			#elif m.count("SV_DEL") > 0:
				#print m
			#	deletions.append(m.split("_")[2])
			elif m.count("SNV") > 0:
				pos = m.split("_")[2]
				start = pos.split(",")[0]
				SNV.append(start)
		status = 4		
		for d in deletions:
			d1 = int(d.split(",")[0])
			d2 = int(d.split(",")[1])
			for s in SNV:
				s1 = int(s)
				if s1 > d1 and s1 < d2:
					status = 5
		return status

filein = open("/ibios/co02/chen/pcawg_2018/result/matrix/coding/select/mutation.csv","r")
fout = open("/abi/data/chen/project_august/mutation2.csv","w")
for line in filein:
	if line[0]=="g":
		line = line.strip("\n")
		#elements = line.split("\t")
		#coh_size = len(elements)-4
		print >>fout,line
		print line
	else:
		line = line.strip("\n")
		elements = line.split("\t")
		cats = []
		for element in elements[1:]:
			dhs = dh_detect(element)
			cats.append(str(dhs))
		print len(cats)
		line2 = "\t".join(cats)		
		print >>fout,elements[0]+"\t"+line2
fout.close()
filein.close()		
