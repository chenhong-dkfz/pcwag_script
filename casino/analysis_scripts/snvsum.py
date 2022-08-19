import sys

filein = sys.argv[1]
fileout = sys.argv[2]

ffi = open(filein,"r")
ffo = open(fileout,"w")
total_length = 0


for line in ffi:
        if line[0]!="#":
	  line = line.strip("\n")
	  elements = line.split("\t")
	  start = int(elements[1])
	  end = int(elements[2])
	  length  = 1	
	  total_length += length
print >>ffo,str(total_length)
