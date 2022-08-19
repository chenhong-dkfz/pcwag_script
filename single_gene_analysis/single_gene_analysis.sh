gene="TP53" # the gene which is interesting
code_dir="/ibios/co02/chen/pcawg_2018/code" # code directory
out_dir="/home/hongc/test" # directory where output PDF file is saved

python ${code_dir}/data_preparation.py ${gene}
/ibios/tbi_cluster/13.1/x86_64/bin/Rscript-3.3.0  ${code_dir}/gene_analysis.R ${gene} ${out_dir}
