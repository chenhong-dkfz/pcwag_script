#for uniuqe the positions and integrate information, separated with ';'
import sys

input_file = open(sys.argv[1],'r') #input file name with argv
output_file = open(sys.argv[2],'w') #output file name with argv
option = sys.argv[3]

#format should be "chr	pos1	pos2	information"

if option == "1":

  chromo0 = '#chromosome'
  pos1 = 'start'
  pos2 = 'end'
  information = 'information'

  for line in input_file:
    line = line.strip('\n')
    elements = line.split('\t')
    if elements[0] != chromo0: #if it is a different chromosome
      if chromo0 != '#chromosome': # if it is not the first description row
	print >>output_file,str(chromo0)+'\t'+str(pos1)+'\t'+str(pos2)+'\t'+information #print current information
      chromo0 = elements[0] #then get new info in
      pos1 = elements[1]
      pos2 = elements[2]
      gene = elements[3]
      if elements[7]!='.': 
	information = elements[7]
      else:
	information = '.'
    elif elements[1] != pos1 or elements[2] != pos2: #if it is a same chromosome, but a new start position
      print >>output_file,str(chromo0)+'\t'+str(pos1)+'\t'+str(pos2)+'\t'+str(gene)+"\t"+information #output what we have
    #print >>output_file,information
      pos1 = elements[1] #then get new info
      pos2 = elements[2]
      gene = elements[3]
      if elements[7]!='.':
	information = elements[7]
      else:
	information = '.'
    elif elements[1]==pos1 and elements[2]==pos2: #if it is not a new region
      
      if elements[7]!='.' and information!='.': #add info or create info
	information = information + ';' + elements[7]
      elif elements[7]!='.' and information=='.':
	information = elements[7]
      else:
	pass
  #else:
  #  print 'error occurs in '+chromo0+' '+' '+elements[1]+' '+pos1+' '+elements[2]+' '+pos2
    
    
  print >>output_file,str(chromo0)+'\t'+str(pos1)+'\t'+str(pos2)+'\t'+str(gene)+"\t"+information 
#print >>output_file,information

if option =="2":
  chromo0 = '#chromosome'
  pos1 = 'start'
  pos2 = 'end'
  gene0 = 'gene'
  score0 = '0'
  count = 1
  for line in input_file:
    line = line.strip("\n")
    elements = line.split("\t")
    if elements[4] =="." or elements[4] == 'score':
      elements[4] = 0
    if elements[0]!= chromo0 or elements[1]!= pos1 or elements[2]!=pos2:
      score0 = str(float(score0)/count)
      if chromo0 != '#chromosome':
	print >>output_file,str(chromo0)+'\t'+str(pos1)+'\t'+str(pos2)+'\t'+str(gene0)+"\t"+str(score0)
      chromo0 = elements[0]
      pos1 = elements[1]
      pos2 = elements[2]
      gene0 = elements[3]
      score0 = elements[4]
      count = 1
    else:
      score0 = str(float(score0)+float(elements[4]))
      count += 1
  score0 = str(float(score0)/count)
  print >>output_file,str(chromo0)+'\t'+str(pos1)+'\t'+str(pos2)+'\t'+str(gene0)+"\t"+str(score0)