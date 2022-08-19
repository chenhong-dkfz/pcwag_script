#for uniuqe the positions and integrate information, separated with ';'
import sys

input_file = open(sys.argv[1],'r') #input file name with argv
output_file = open(sys.argv[2],'w') #output file name with argv

#input_file = open("/ibios/co02/chen/casino/promoter/intermediate/PRAD-US/252a8545-fdc7-42d5-9041-c459b44095b9/snv_tmp_rsv.bed","r")
#output_file = open("/home/hongc/test/sample_output.txt","w")



data = {}

for line in input_file:
  line = line.strip("\n")
  elements = line.split("\t")
  chromosome = elements[0]
  pos1 = elements[1]
  pos2 = elements[2]
  gene = elements[3]
  information = elements[7]
  position = pos1+"\t"+pos2+"\t"+gene
  #if information != ".":
  infos = information.split(";")
  if position in data.keys():
    data[position] += infos
  else:
    data[position] = infos
    
for item in sorted(data.keys()):
  for val in data[item]:
    if val == ".":
      data[item].remove(val)
  print >>output_file,item+"\t"+str(data[item])