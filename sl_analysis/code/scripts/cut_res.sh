path=/home/hongc/test/result_subcohort2/all
path2=/home/hongc/test/result_subcohort22/all
mkdir -p ${path2}

for f in ${path}/*.txt ; do 
  echo ${f}
  #cut -f 1,4,5,7,8,9,10,11,12 
done