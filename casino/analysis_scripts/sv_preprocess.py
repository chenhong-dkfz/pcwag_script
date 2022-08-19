import sys
import os

coh = sys.argv[1]
bed_path = sys.argv[2]
new_bed_path = bed_path+"/"+coh

if not os.path.isdir(new_bed_path):
  os.mkdir(new_bed_path)


match_list = open("/home/hongc/CaSINo/coh_pid/transfer.txt","r")

match_dict = {}

for line in match_list:
  line = line.strip("\n")
  ee = line.split("\t")
  match_dict[ee[0]]=ee[1]

match_list.close()

pids = open("/home/hongc/CaSINo/coh_pid/"+coh+".txt","r")

for line in pids:
  line = line.strip("\n")
  new_id = match_dict[line]
  
  #print new_id
  bed_name = bed_path+"/"+new_id+".sv.bed"
  if os.path.isfile(bed_name):
    os.system("cp "+bed_name+" "+new_bed_path)

pids.close()


