
def data_redo(idx,cohort):
	string = datas[idx]
	string = string.split("\t")
	selected = [i for i, x in enumerate(coh_line) if x == cohort]
	string = [string[i] for i in selected]
	string = "\t".join(string)
	return(string)


idx = 1
datas = ["1\t2\t3\t4\t5\t6","7\t8\t9\t10\t11\t12"]
coh_line = ["a","a","b","a","a","b"]

c = data_redo(idx,"a")
print c

