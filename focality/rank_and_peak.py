chrs = []
starts = []
ends = []
genes = []

fin = open("/abi/data/chen/co02/pcawg_2018/result/matrix/coding/LGG-US2.bed","r")
for line in fin:
	if line[0]!="#":
		line = line.strip("\n")
		el = line.split("\t")
		ch = el[0]
		start = el[1]
		end = el[2]
		gene = el[3]
		chrs.append(ch)
		starts.append(start)
		ends.append(end)
		genes.append(gene)





# return a list that if the gene is at the edge, the index is 1=start, 2=end, otherwise is 0
end_index = [0]*len(chrs)
for i in range(0,len(chrs)):
	if i == 0 or chrs[i]!=chrs[i-1]:
		end_index[i] = 1
	elif i == (len(chrs)-1):
		end_index[i] = 2
	elif chrs[i]!=chrs[i+1]:
		end_index[i] = 2

def get_nb(x): #default 5 in both flanks
	min_index = 0
	max_index = len(chrs)
	if x-5 < min_index:
		list1 = range(min_index,x+6)
	elif x+5 >= max_index:
		list1 = range(x-5,max_index)
	else:
		list0 = [x-5,x-4,x-3,x-2,x-1,x,x+1,x+2,x+3,x+4,x+5]
		chr_list = [chrs[i] for i in list0]
		match_list = [i for i, p in enumerate(chr_list) if p == chrs[x]]
		list1 = [list0[i] for i in match_list]
	return list1	

def get_range(x,window):
	start0 = starts[x]
	end0 = ends[x]
	
	list1 = []
	list1.append(x)
	if x!=0 and x!=len(chrs)-1:
		i = 1
		while int(start0)-int(ends[x-i])<=window and int(start0)-int(ends[x-i])>0 and x-i>=0:
			list1.append(x-i)
			i += 1
		i = 1
		while int(starts[x+i])-int(end0)<=window and int(starts[x+i])-int(end0)>0 and x+i<len(chrs):
			list1.append(x+i)
			i += 1
		list1 = list(set(list1))	
		list1 = sorted(list1)
	elif x==0:
		i = 1
		while int(starts[x+i])-int(end0)<=window and int(starts[x+i])-int(end0)>0 and x+i<len(chrs):
			list1.append(x+i)
			i += 1
		list1 = list(set(list1))	
		list1 = sorted(list1)
	elif x==len(chrs)-1:
		i = 1
		while int(start0)-int(ends[x-i])<=window and int(start0)-int(ends[x-i])>0 and x-i>=0:
			list1.append(x-i)
			i += 1
		list1 = list(set(list1))	
		list1 = sorted(list1)	
	return(list1)

# the rank is defined as 0,1,2,....
def get_rank(list_k,x):
	list_k = sorted(list_k,reverse=1)
	rank = list_k.index(x)
	return rank


def peak_score(x,sc):
	if end_index[x]==0:
		peak_score = sc[x]*2-sc[x-1]-sc[x+1]
	elif end_index[x]==1:
		peak_score = sc[x]*2-sc[x+1]*2
	elif end_index[x]==2:
		peak_score = sc[x]*2-sc[x-1]*2
	return peak_score


# read fscores
fs_6 = []
fin_fs = open("/abi/data/chen/newhome/promoter_focallity/data/all_fscore_1e6.txt","r")
for line in fin_fs:
	sc = line.strip("\n")
	sc = float(sc)
	fs_6.append(sc)


fs_7 = []
fin_fs = open("/abi/data/chen/newhome/promoter_focallity/data/all_fscore_1e7.txt","r")
for line in fin_fs:
	sc = line.strip("\n")
	sc = float(sc)
	fs_7.append(sc)






# return a list of neighbours of gene, two ways, 5 genes in each flank which contain 11 genes in total; 1e7 in each flank which contain several genes
nb_list_1 = []
range_list_1 = []

for i in range(0,len(chrs)):
	nb_list_1.append(get_nb(i))
	range_list_1.append(get_range(i,10000000))


print(len(nb_list_1))
print(len(range_list_1))
print(nb_list_1[0])
for i in range(1,10):
	print(range_list_1[0])		

ftest = open("/abi/data/chen/newhome/promoter_focallity/data/intermediate/int2.tsv","w")
fnb = open("/abi/data/chen/newhome/promoter_focallity/data/intermediate/nb_info2.tsv","w")

for i in range(0,len(chrs)):
	print i
	nb_fs6_list = [fs_6[j] for j in nb_list_1[i]]
	nb_fs6_rank = get_rank(nb_fs6_list,fs_6[i]) 
	nb_fs7_list = [fs_7[j] for j in nb_list_1[i]]
	nb_fs7_rank = get_rank(nb_fs7_list,fs_7[i])
	range_fs6_list = [fs_6[j] for j in range_list_1[i]]
	range_fs6_rank = get_rank(range_fs6_list,fs_6[i])
	range_fs7_list = [fs_7[j] for j in range_list_1[i]]
	range_fs7_rank = get_rank(range_fs7_list,fs_7[i])
	peak_nb_fs6 = peak_score(i,fs_6)
	peak_nb_fs7 = peak_score(i,fs_7)
	print >>fnb,str(chrs[i])+"\t"+str(starts[i])+"\t"+str(ends[i])+"\t"+str(nb_list_1[i])+"\t"+str(range_list_1[i])+"\t"+str(i)
	print >>ftest,str(fs_6[i])+"\t"+str(fs_7[i])+"\t"+str(nb_fs6_rank)+"\t"+str(nb_fs7_rank)+"\t"+str(range_fs6_rank)+"\t"+str(range_fs7_rank)+"\t"+str(peak_nb_fs6)+"\t"+str(peak_nb_fs7)

ftest.close()
fnb.close()