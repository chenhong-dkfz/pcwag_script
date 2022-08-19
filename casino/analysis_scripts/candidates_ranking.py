import sys
import operator

#input_file = open("/ibios/co02/chen/casino/results/analysis/qqplot_PRAD-US.txt","r") #input combined functional elements with scores
#output_file = open("","w") #output 
result_dir = sys.argv[1]
coh = sys.argv[2]
filein = result_dir+"/analysis/candidates_"+coh+".txt"
fileout = result_dir+"/analysis/selected_candidates_"+coh+".txt"
input_file = open(filein,"r") #input combined functional elements with scores
output_file = open(fileout,"w") #output 


candidates = {}
snv_score = {}
var_score = {}

for line in input_file:
  if line[0]=="c" or line[0] == "#":
    header = line.strip()
  else:
    line = line.strip()
    elements = line.split("\t")
    region = "\t".join(elements[0:3])+"\t"+elements[4]
    
    
    if region not in candidates.keys():
      candidates[region] = 1
    else:
      candidates[region] += 1
     
    if region not in snv_score.keys():
      snv_score[region] = float(elements[3])
    else:
      snv_score[region] += float(elements[3])
      
    if region not in var_score.keys():
      var_score[region] = float(elements[6])
    else:
      var_score[region] += float(elements[6])
   
#print candidates   
sorted_candidates = sorted(candidates.items(),key=operator.itemgetter(1),reverse=True)
for s in sorted_candidates:
  print >>output_file,s[0]+"\t"+str(s[1])+"\t"+str(snv_score[s[0]]/s[1])+"\t"+str(var_score[s[0]]/s[1])