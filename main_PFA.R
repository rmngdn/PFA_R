#MAIN PFA
#DONT FORGET THE REFERENCES 
#AND TO PACKAGE IT IN THE END

source("./PFA_algorithms.R")



## Set input and output data filename:

# for example:
# input:
gene_fileName <- "data_gene_expression.csv"
methy_fileName <- "data_methy_expression.csv"
miRNA_fileName <- "data_mirna_expression.csv"
fileName_list <- list(gene_fileName, 
                      methy_fileName, 
                      miRNA_fileName)

# output:
output_fileName <- "global_sample_spectrum.csv"

## Import the data:

# for example
data_list <- lapply(fileName_list, read.csv)

## Capture the local sample-spectrum for each biological data type

