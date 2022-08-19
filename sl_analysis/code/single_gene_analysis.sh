gene="RBL2" # the gene which is interesting
code_dir="/home/hongc/PhD/TumorPrint/single_gene/code" # the code directory contains a python, a R as well as a shell script
out_dir="/home/hongc/test/result_080818" # directory where output PDF file saves

mkdir -p ${out_dir}
python ${code_dir}/data_preparation.py ${gene}
/ibios/tbi_cluster/13.1/x86_64/bin/Rscript-3.3.0  ${code_dir}/gene_analysis.R ${gene} ${out_dir}
