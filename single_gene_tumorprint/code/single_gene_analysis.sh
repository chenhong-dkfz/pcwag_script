gene="USMG5" # the gene which is interesting
code_dir="/abi/data/chen/newhome/Tumorprint/single_gene/code" # the code directory contains a python, a R as well as a shell script
out_dir="/abi/data/chen/newhome/Tumorprint/single_gene/result" # directory where output PDF file saves

mkdir -p ${out_dir}
python ${code_dir}/data_preparation.py ${gene}
/ibios/tbi_cluster/13.1/x86_64/bin/Rscript-3.3.0  ${code_dir}/gene_analysis.R ${gene} ${out_dir}
