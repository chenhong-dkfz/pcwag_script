import os


for file in os.listdir("/ibios/co02/chen/pcawg_2018/snv_exon"):
    if file.endswith(".txt"):
        fin = open(os.path.join("/ibios/co02/chen/pcawg_2018/snv_exon", file),"r")
        fout = open(os.path.join("/ibios/co02/chen/pcawg_2018/snv_exon2", file),"w")
        for line in fin:
        	line = line.strip("\n")
        	el = line.split("\t")
        	snv = el[3]+"_"+el[1]+","+el[2]
        	print >>fout, el[0]+"\t"+el[1]+"\t"+el[2]+"\t"+snv
        fin.close()
        fout.close()
