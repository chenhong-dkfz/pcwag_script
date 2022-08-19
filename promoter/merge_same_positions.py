import sys

input_file = open(sys.argv[1],'r') #input file name with argv
output_file = open(sys.argv[2],'w') #output file name with argv
#option = sys.argv[3]

chromo0 = '#chromosome'
pos1 = 'start'
pos2 = 'end'
gene = 'gene'
information = 'information'

for line in input_file:
	line = line.strip("\n")

	if len(line) < 2:
		pass
	else:
		elements = line.split('\t')
		#print elements
		if elements[0] != chromo0:
			if chromo0 != '#chromosome':
				#print str(chromo0)+'\t'+str(pos1)+'\t'+str(pos2)+'\t'+str(gene)+'\t'+information
				print >>output_file,str(chromo0)+'\t'+str(pos1)+'\t'+str(pos2)+'\t'+str(gene)+'\t'+information
			chromo0 = elements[0]
			pos1 = elements[1]
			pos2 = elements[2]
			gene = elements[3]
			if elements[7]!='.':
				information = elements[7]
			else:
				information = '.'

		elif elements[1] != pos1 or elements[2] != pos2 or elements[3]!=gene:
			print >>output_file,str(chromo0)+'\t'+str(pos1)+'\t'+str(pos2)+'\t'+str(gene)+"\t"+information
			#print str(chromo0)+'\t'+str(pos1)+'\t'+str(pos2)+'\t'+str(gene)+"\t"+information
			pos1 = elements[1]
			pos2 = elements[2]
			gene = elements[3]
			if elements[7]!='.':
				information = elements[7]
			else:
				information = "."
		elif elements[1]==pos1 and elements[2]==pos2 and elements[3]==gene:
			if elements[7]!='.' and information!='.':
				information = information + ';' + elements[7]
			elif elements[7]!='.' and information=='.':
				information = elements[7]
			else:
				pass
if chromo0 != '#chromosome':
	#print str(chromo0)+'\t'+str(pos1)+'\t'+str(pos2)+'\t'+str(gene)+"\t"+information
	print >>output_file,str(chromo0)+'\t'+str(pos1)+'\t'+str(pos2)+'\t'+str(gene)+"\t"+information		

input_file.close()
output_file.close()
