#for new version output bed file

import sys

input_file = open(sys.argv[1],'r') #input file name with argv
output_file = open(sys.argv[2],'w') #output file name with argv

current_chromo = "CC"
current_1 = "0"
current_2 = "1"
current_gene = ""
count = 0
count2 = 0
count3 = 0

for line in input_file:
  count3 += 1
  line = line.strip("\n")
  elements = line.split("\t")
  chromosome = elements[0]
  pos1 = elements[1]
  pos2 = elements[2]
  gene = elements[3]
  if current_chromo == chromosome and current_1 == pos1 and current_2 == pos2:
    current_gene = current_gene+";"+gene
    count += 1
  else:
    if current_chromo != "CC":
      print >> output_file,current_chromo+"\t"+current_1+"\t"+current_2+"\t"+current_gene
      count2 += 1
    current_chromo = chromosome
    current_1 = pos1
    current_2 = pos2
    current_gene = gene
      
      

print >> output_file,current_chromo+"\t"+pos1+"\t"+pos2+"\t"+gene    
count2 += 1

print count
print "count2 "+str(count2)
print "count3 "+str(count3)