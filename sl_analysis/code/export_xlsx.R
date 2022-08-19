# library("tidyverse")
library("cli")
library("writexl")

setwd("/home/hongc/test/result_subcohort2")
#setwd("/abi/data/chen/project_august/mutation25/cohortwise")


FileName_import <- NULL
#FileName_import <- "result/pten.rdata"
FileName_import <- "result/alc1.rdata"

FileNameExport <- NULL
FileNameExport <- sub(x=FileName_import, pattern="\\.rdata$", replacement=".xlsx")

stopifnot(file.exists(FileName_import))


legend_descriptor <- NULL
legend_descriptor <- data.frame(
    label=c(
        ## selected subcohort
        "dataset",
        "histology",
        "gene_of_interest",
        ## the names of the tables (xlsx sheet names corresponding to the 25 comparisons), e.g.:
        "mut_hiCo",
        "mut_medCo",
        "mut_lowCo",
        "amp_hi",
        "amp_norm",
        ## table fields
        "chr",
        "gene_name",
        "cgc",
        "p",
        "OR",
        "sample_size",
#         "PTEN_mutations",
        "CHD1L_mutations",
        "gene_mutation",
        "overlaps"
    ),
    description=c(
        ## selected subcohort
        "BRCA-EU, BRCA-US, ...",
        "breast, breast, ...",
#         "PTEN",
        "CHD1L",
        ## the names of the tables (xlsx sheet names corresponding to the 25 comparisons), e.g.:
        "bi-allelic mutations only",
        "bi-allelic + potential bi-allelic mutations",
        "bi-allelic + potential bi-allelic + mono-allelic mutations + ",
        "copy number gains >8",
        "copy number gains >4",
        ## table fields
        "chromosome the gene is annotated to",
        "official gene symbol by HUGO nomenclature",
        "cancer gene census annotation",
        "p-value of hypergeometric test",
        "effect strength (odds ratio)",
        "total number of samples in the selected subcohort(s)",
#         "count of samples with PTEN mutation",
        "count of samples with CHD1L mutation",
        "count of samples with mutation in other gene",
        "count of samples with mutation in both genes"
        ),
    stringsAsFactors=FALSE
)

datalist <- NULL
load(FileName_import)

names_vector <- NULL
names_vector <- c("mut_lowCo", "amp_norm", "mut_medCo", "amp_hi", "mut_hiCo")
#names_vector <- c("mut_hiCo", "mut_medCo", "mut_lowCo", "amp_hi", "amp_norm")
# names_vector <- as.vector(outer(X=names_vector, Y=names_vector, FUN=paste, sep=" - ")) # by column
names_vector <- as.vector(t(outer(X=names_vector, Y=names_vector, FUN=paste, sep=" - "))) # by row

stopifnot(length(datalist) == length(names_vector))
names(datalist) <- names_vector

datalist[["legend"]] <- legend_descriptor

## re-order legend to beginning
datalist <- datalist[c("legend", names_vector)]



cli::cat_line("Writing ", length(datalist), " tables to file '", FileNameExport, "'.")
write_xlsx(x=datalist, path="/home/hongc/test/alc1.xlsx", col_names=TRUE)


