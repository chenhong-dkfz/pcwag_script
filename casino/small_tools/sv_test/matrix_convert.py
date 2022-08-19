import re
import sys

filein = open(sys.argv[1],"r")
fout = open(sys.argv[2],"w")

for line in filein:
  if line[0]=="#":
    line = line.strip("\n")
    print >>fout,line,
  else:
    line = line.strip("\n")
    elements = line.split("\t")
    chromo = elements[0]
    start = elements[1]
    end = elements[2]
    gene = elements[3]
    print >>fout,chromo+"\t"+start+"\t"+end+"\t"+gene+"\t",
    for element in elements[4:]:
      #outsv = []
      svs = re.findall("DEL_\d+",element)
      min_len = "default"
      for sv in svs:
        m = re.match("DEL_(\d+)",sv)
        if m:
          if min_len == "default":
            min_len = int(m.groups()[0])
          #min_sv = sv
          elif min_len > int(m.groups()[0]):
            #min_sv = sv
            min_len = int(m.groups()[0])
      if min_len == "default":
        min_len = 0
      element = str(min_len)
      print >>fout,element+"\t",
  print >>fout,""
