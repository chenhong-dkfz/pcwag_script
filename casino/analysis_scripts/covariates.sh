cut -f 1,2,3 ${mutation_file} > ${coordinate}
bedtools intersect ${coordinate}