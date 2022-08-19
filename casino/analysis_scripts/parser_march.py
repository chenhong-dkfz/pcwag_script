


#snvs_ICGC_MB38_somatic_snvs_conf_8_to_10.vcf
#print >>fo,chrom+"\t"+pos+"\t"+end+"\t"+"SNV_"+ref+"_"+alt+"_"+str(af)+"_"+pos
import os
import os.path
import re
import sys

pid = sys.argv[1]
coh = sys.argv[2]
output_path = sys.argv[3]

#vcf_path = '/icgc/lsdf/mb/analysis/medullo/reanno_ivo/pancanworkflow/results_per_pid_all/'+pid+'/mpileup/'
#vcf_file = vcf_path + "snvs_"+pid+"_somatic_snvs_conf_8_to_10.vcf"

#vcf_path = '/ibios/co02/chen/data_march_2016/variants_per_cohorts/'+coh+'/'
vcf_path = '/ibios/co02/chen/MB/MB_vcfs/'
vcf_file = vcf_path + "snvs_"+pid+"_somatic_snvs_conf_8_to_10.vcf"


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


fin = open(vcf_file,"r")
fo = open(parser_path+"snv.bed","w")


for line in fin:
  if line[0]!="#":
    line = line.strip("\n")
    elements = line.split("\t")
    chrom = elements[0]
    pos = elements[1]
    end = str(int(pos)+1)
    ref = elements[3]
    alt = elements[4]
    #fil = elements[6]
    #ctx = elements[10]
    info = elements[7]
    #s = re.match(r'.*?medianVAF=(.*)',info) 
    #vaf= s.groups()[0]
    

    #if fil!="PASS":
    #  pass
    #else:
    #  end = str(int(pos)+1)
    #  n = re.match(r'\w+?(\w),(\w)\w+',ctx)
    #  pos_1 = n.group(1)
    #  pos_2 = n.group(2)
  
      #catogory = 0
      #if (ref=="A" and alt=="G") or (ref=="T" and alt=="C"):
#	catogory = 5 #AT transition
#      elif (ref=="A" and (alt=="C" or alt=="T")) or (ref=="T" and (alt=="A" or alt=="G")):
#	catogory = 6 #AT transversion
#      elif (ref=="C" and alt=="T") or (ref=="G" and alt=="A"):
#	if (ref=="C" and pos_2=="G") or (ref=="G" and pos_1=="C"):
#	  catogory = 1 #CpG transition
#	else:
#	  catogory = 3 #CG transition
 #     elif (ref=="C" and (alt=="A" or alt=="G")) or (ref=="G" and (alt=="C" or alt=="T")):
#	if (ref=="C" and pos_2=="G") or (ref=="G" and pos_1=="C"):
#	  catogory = 2 #CpG transversion
#	else:
#	  catogory = 4 #CG transversion
      
#    print >>fo,chrom+"\t"+pos+"\t"+end+"\t"+"SNV_"+ref+"_"+alt+"_"+vaf+"_"+pos
    print >>fo,chrom+"\t"+pos+"\t"+end+"\t"+"SNV_"+ref+"_"+alt+"_"+pos
