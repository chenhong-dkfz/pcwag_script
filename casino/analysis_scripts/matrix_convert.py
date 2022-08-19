#this script is used to convert original matrix into required forms

import re
import sys

with_cnv = sys.argv[1]
fi = sys.argv[2]
fo1 = sys.argv[3]
fo2 = sys.argv[4]

filein = open(fi,"r")
fileout1 = open(fo1,"w")


#filein = open("/home/hongc/test.bed","r")
#fileout = open("/home/hongc/test_snv.bed","w")


#options are: 1.snv only;2.snv only frequency;3.cnv only; 4.cnv only minimail length
option = 2
cnv_option = 4
#additional option: a.cnv all, d.cnv deletion, r.cnv replication, l,cnv LOH

additional_option = "d"

#this part is for snv only

if option == 1 or option ==2:
	for line in filein:
		if line[0]=="#":
			line = line.strip("\n")
			print >>fileout1,line
		else:
			line = line.strip("\n")
			elements = line.split("\t")
			chromo = elements[0]
			start = elements[1]
			end = elements[2]
			gene = elements[3]
			print >>fileout1,chromo+"\t"+start+"\t"+end+"\t"+gene+"\t",
			for element in elements[4:]:
				if option == 1:
					element = ";".join(re.findall("SNV_\w+_\w+_\d+\.\d+_\d+",element))
					if len(element)==0:
						element = "."
					print >>fileout1,element+"\t",
				if option == 2:
					element = element.count("SNV")
					element = str(element)
					print >>fileout1,element+"\t",
			print >>fileout1,""

#this part is for cnv only
fileout1.close()
filein = open(fi,"r")

if with_cnv == "yes":
	fileout2 = open(fo2,"w")
	if cnv_option == 3 or cnv_option == 4:
		for line in filein:
			if line[0]=="#":
				line = line.strip("\n")
				print >>fileout2,line
			else:		
				line = line.strip("\n")
				elements = line.split("\t")
				chromo = elements[0]
				start = elements[1]
				end = elements[2]
				gene = elements[3]
				print >>fileout2,chromo+"\t"+start+"\t"+end+"\t"+gene+"\t",
				for element in elements[4:]:
					if cnv_option == 3:  
						outcnv = [] 
						#element = ";".join(re.findall("CNV_\d+M_\d+_\d+_\d+",element))
						allcnv = re.findall("CNV_\d+M_\d+_\d+_\d+",element)
						for onecnv in allcnv:
							m = re.match("CNV_\d+M_(\d+)_(\d+)_(\d+)",onecnv)
							if additional_option == "d": #this option is for cnv deletions
								if (m.groups()[1]=="1" or m.groups()[2]=="0"):
									outcnv.append(onecnv)
							#then if all,replcation and loss of heretogolous
							elif additional_option == "r":
								if (m.groups()[1]>"2" or m.groups()[2]!="0"):
									outcnv.append(onecnv)
							elif additional_option == "l":
								if (m.groups()[1]>"1" or m.groups()[2]=="0"):
									outcnv.append(onecnv)
							elif additional_option == "a":
								outcnv.append(onecnv)
							else:
								print "invalid options!"


						element = ";".join(outcnv)			
						if len(element)==0:
							element = "."
						print >>fileout2,element+"\t",
					if cnv_option == 4:
						#m = re.match("CNV_\d+M_\d+_\d+_\d+",element)
						outcnv = [] 
						cnvs = re.findall("CNV_\d+M_\d+_\d+_\d+",element)
						min_len = "default"
						for onecnv in cnvs:
							m = re.match("CNV_\d+M_(\d+)_(\d+)_(\d+)",onecnv)
							if additional_option == "d": #this option is for cnv deletions
								if (m.groups()[1]=="1" or m.groups()[2]=="0"):
									outcnv.append(onecnv)
									#print outcnv
							#then if all,replcation and loss of heretogolous
							elif additional_option == "r":
								if (m.groups()[1]>"2" or m.groups()[2]!="0"):
									outcnv.append(onecnv)
							elif additional_option == "l":
								if (m.groups()[1]>"1" or m.groups()[2]=="0"):
									outcnv.append(onecnv)
							elif additional_option == "a":
								outcnv.append(onecnv)
							else:
								print "invalid options!"

						for cnv in outcnv:
							m = re.match("CNV_\d+M_(\d+)_\d+_\d+",cnv)
							if min_len == "default":
								min_len = int(m.groups()[0])
								min_cnv = cnv
							elif min_len > int(m.groups()[0]):
								min_cnv = cnv
								min_len = int(m.groups()[0])
						if min_len == "default":
							min_len = 0		
						element = str(min_len)
						print >>fileout2,element+"\t",
				print >>fileout2,""
	fileout1.close()			



