gene="USMG5" # the gene which is interesting
code_dir="/omics/groups/OE0436/data/chen/home/single_gene_tumorprint" # the code directory contains a python, a R as well as a shell script
out_dir="/omics/groups/OE0436/data/chen/home/test/result_08042022" # directory where output PDF file saves

mkdir -p ${out_dir}
python ${code_dir}/data_preparation.py ${gene}
module load r/3.5.1
#/ibios/tbi_cluster/13.1/x86_64/R/R-3.5.0/bin/Rscript ${code_dir}/gene_analysis.R ${gene} ${out_dir}
Rscript ${code_dir}/gene_analysis.R ${gene} ${out_dir}
