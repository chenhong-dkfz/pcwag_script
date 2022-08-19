import sys
import re
import os
import gzip


file_path = "/abi/data/lars/pcawg/data_march2016/sv_merge_draft_April2016/"
pid_id = sys.argv[1]
file_suffix_1 = ".pcawg6_merge_1-0-0.160419.somatic.sv.vcf.gz"
file_suffix_2 = ".pcawg6_merge_1-0-0.160420.somatic.sv.vcf.gz"


file_name_1 = file_path + pid_id + file_suffix_1
file_name_2 = file_path + pid_id + file_suffix_2

if os.path.isfile(file_name_1):
  file_name = file_name_1
elif os.path.isfile(file_name_2):
  file_name = file_name_2
else:
  print "NO DATA FOR"+pid_id
  

#output_path = "/home/hongc/CaSINo/code/small_tools/sv_test/temp/"
#output_suffix = ".sv_parsed.bed"
#output_name = output_path + pid_id + output_suffix

output_name = sys.argv[2]


pattern1 = re.compile(r"(\w*)\](\d+):(\d+)\]")
pattern2 = re.compile(r"\](\d+):(\d+)\](\w*)")
pattern3 = re.compile(r"(\w*)\[(\d+):(\d+)\[")
pattern4 = re.compile(r"\[(\d+):(\d+)\[(\w*)")
#pattern0 = re.compile(r".*SVCLASS=(\w+);.*SVLEN=(\d+).*")
pattern0 = re.compile(r".*SVCLASS=(\w+);.*")
file_in = gzip.open(file_name,"r")
file_out = open(output_name,"w")

#select_type = "all"
select_type = "DEL"
#select_type = "DUP"


print >>file_out,"#chromo"+"\t"+"start"+"\t"+"end"+"\t"+"sv_class"+"\t"+"chromo_vali"+"\t"+"pid_id"+"\t"+"next_chrom"


for line in file_in:
  if line[0]!="#":
    line = line.strip("\n")
    elements = line.split("\t")
    chromo = elements[0]
    start = elements[1]
    alt = elements[4]
    alt = re.sub('chr','',alt)

    info = elements[7]
    if pattern1.match(alt):
      m = pattern1.match(alt)
      t = m.groups()[0]
      next_chrom = m.groups()[1]
      end = m.groups()[2]
    elif pattern2.match(alt):
      m = pattern2.match(alt)
      t = m.groups()[0]
      next_chrom = m.groups()[1]
      end = m.groups()[2]
    elif pattern3.match(alt):
      m = pattern3.match(alt)
      t = m.groups()[0]
      next_chrom = m.groups()[1]
      end = m.groups()[2]
    elif pattern4.match(alt):
      m = pattern4.match(alt)
      t = m.groups()[0]
      next_chrom = m.groups()[1]
      end = m.groups()[2]      
    else:
      t = "error"
      next_chrom = "None"
      end = "123456789"
    try:
      m = pattern0.match(info)
      typ = m.groups()[0]
     # length = m.groups()[1]
    except:
      typ = "UNKNOWN0"
     # length = "0"
    #if (int(end)-int(start))!=int(length):
    #  length_vali = "F"
    #  print int(end)
    #  print int(start)
    #  print int(length)
    #  print "###"
    #else:
    #  length_vali = "T"
    if chromo == next_chrom:
      chromo_vali = "T"
    else:
      chromo_vali = "F"
    if chromo == next_chrom:
      length = str(int(end)-int(start))
      if int(start)>int(end):
        mid = start
        start = end
        end = mid
        length = str(int(end)-int(start))
    else:
        length = "0"
    if chromo == next_chrom:
      if typ == select_type:
        print >>file_out,chromo+"\t"+start+"\t"+end+"\t"+typ+"_"+length+"\t"+chromo_vali+"\t"+pid_id+"\t"+next_chrom       

  
