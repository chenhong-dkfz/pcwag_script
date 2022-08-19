path="/home/hongc/test/result_subcohort3"
muts="mut nam pam pa ra"
for muti in ${muts};do
  for mutj in ${muts};do
    mut=${muti}_${mutj}
    echo ${mut}
    mkdir -p ${path}/${mut}
    cp ${path}/*${mut}.txt ${path}/${mut}
    #echo ${path}/*${mut}.txt
    
  done
done
  
