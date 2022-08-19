import gzip
import re
import sys
import os
import os.path
import pybedtools

pid = sys.argv[1]
coh = sys.argv[2]
output_path = sys.argv[3]
#output_path = sys.argv[3]
#vcf_gz_path = sys.argv[4]

vcf_path = '/icgc/pcawg/project/results/variants/sanger/data/'+pid+'/'+pid+'/'

vcffiles = [f for f in os.listdir(vcf_path) if os.path.isfile(os.path.join(vcf_path,f))]
#output_path = '/ibios/co02/chen/enhancer/'
#output_path = '/ibios/co02/chen/promoter/' #alternative for promoter


#create directory if necessary

if not os.path.exists(output_path+"/"+coh):
  try: 
    os.mkdir(output_path+"/"+coh)
  except:
    pass
if not os.path.exists(output_path+"/"+coh+'/'+pid):
  try: 
    os.mkdir(output_path+"/"+coh+'/'+pid)
  except:
    pass
parser_path = output_path+"/"+coh+'/'+pid+'/'

#ee = open('/ibios/co02/chen/tmp/'+pid+'cnv_overlaps.txt',"a") ##here we check if there are overlaps between cnvs

for vcffile in vcffiles:
  try:
    match = re.match(".*?\.somatic\.(.*)\.vcf\.gz$",vcffile)

    if match.groups(1)[0]=="sv":
      sv = os.path.join(vcf_path,vcffile)
    elif match.groups(1)[0]=="cnv":
      cnv = os.path.join(vcf_path,vcffile)
    elif match.groups(1)[0]=="indel":
      indel = os.path.join(vcf_path,vcffile)
    elif match.groups(1)[0]=="snv_mnv":
      snv = os.path.join(vcf_path,vcffile)
  except:
    continue

if not os.path.isfile(parser_path+"cnv_test_pass.bed"):
  table = gzip.open(cnv)
  fo = open(parser_path+"cnv_test_pass.bed","w")
  last_end = 0
  last_chrom = "1"
  for line in table:
      if line[0]!="#":
          elements = line.split("\t")
          chrom = elements[0]
          pos = elements[1]
          
          if chrom!=last_chrom:
	    last_chrom = chrom
	    last_end = 0
          
          #if int(pos)<last_end:
	  #  print>>ee,pid
	    
          info = elements[7]
          m = re.match(r'SVTYPE=(.+?);END=(\d+)',info)
          typ = m.group(1)
          end = m.group(2)
          
          last_end = int(end)
          
          length = str(int(end)-int(pos)) #here without +1
          if int(length)>1000000000:
            level = "1000M"
          elif int(length)>100000000:
            level = "100M"
          elif int(length)>10000000:
            level = "10M"
          elif int(length)>1000000:
            level = "1M"
          elif int(length)>100000:
            level = "100K"
          elif int(length)>10000:
            level = "10K"
          else:
            level = "1K"
          tumor = elements[10]
          l = re.match(r'.*(\d+):(\d+$)',tumor)
          cn = l.group(1)
          cn2 = l.group(2)
          if int(cn)==2 and int(cn2)==1:
	    pass
	  else:
	    print >>fo,chrom+"\t"+pos+"\t"+end+"\t"+typ+"_"+level+"_"+str(length)+"_"+cn+"_"+cn2
  table.close()
  fo.close()

#indel files#
if not os.path.isfile(parser_path+"indel_test_pass.bed"):
  table = gzip.open(indel)
  fo = open(parser_path+"indel_test_pass.bed","w")

  for line in table:
      if line[0]!="#":
          elements = line.split("\t")
          chrom = elements[0]
          pos = elements[1]
          ref = elements[3]
          length = len(ref)
          fil = elements[6]
          if fil!='PASS':
	    pass
	  else:
	    end = str(int(pos) + length) #here without -1
	    print >>fo,chrom+"\t"+pos+"\t"+end+"\t"+"INDEL"
  table.close()
  fo.close()


#snv_mnv files#
if not os.path.isfile(parser_path+"snv_mnv_test_pass.bed"):
  table = gzip.open(snv)
  fo = open(parser_path+"snv_mnv_test_pass.bed","w")

  for line in table:
      if line[0]!="#":
          elements = line.split("\t")
          chrom = elements[0]
          pos = elements[1]
          ref = elements[3]
          alt = elements[4]
          fil = elements[6]
          tumor = elements[10]
          if fil!='PASS':
	    pass
	  else:
	    end = str(int(pos) + len(ref)) #here without -1
          #print tumor
	    n = re.match(r'^\d+\|(\d+):(\d+):(\d+):(\d+):(\d+):(\d+):(\d+):(\d+):(\d+):.*',tumor)
          #print n.group(3)
	    su = 0
	    for i in range(2,10):
		su = su + int(n.group(i))
          #print su
	    if alt=="A":
		sa = int(n.group(2)) + int(n.group(6))
	    elif alt=="C":
		sa = int(n.group(3)) + int(n.group(7))
	    elif alt=="G":
		sa = int(n.group(4)) + int(n.group(8))
	    elif alt=="T":
		sa = int(n.group(5)) + int(n.group(9))
	    af = float(sa)/float(su)
	    af = round(af,2)
	    print >>fo,chrom+"\t"+pos+"\t"+end+"\t"+"SNV_"+ref+"_"+alt+"_"+str(af)+"_"+pos
  table.close()
  fo.close()

#sv_files#
if not os.path.isfile(parser_path+"sv_test_pass.bed"):
  table = gzip.open(sv)
  fo = open(parser_path+"sv_test_pass.bed","w")

  for line in table:
      if line[0]!="#":
          elements = line.split("\t")
          chrom = elements[0]
          pos = elements[1]
          info = elements[7]
          p = re.match(r'SVTYPE=(.+?);.*',info)
          typ = p.group(1)
          end = str(int(pos)+1)
          print >>fo,chrom+"\t"+pos+"\t"+end+"\t"+"SV_"+typ
  table.close()
  fo.close()

#ee.close()
