#MAIN PFA
#DONT FORGET THE REFERENCES 
#AND TO PACKAGE IT IN THE END

pathToYourWD <- "~/Documents/Local_Work/PFA to R"
setwd(pathToYourWD)
source("PFA_algorithms.R")



## Set input and output data filename:

# for example:
# input:
gene_fileName <- "GLIO_Gene_Expression.txt"
methy_fileName <- "GLIO_Methy_Expression.txt"
miRNA_fileName <- "GLIO_Mirna_Expression.txt"
fileName_list <- list(gene_fileName, 
                      methy_fileName, 
                      miRNA_fileName)

# output:
output_fileName <- "global_sample_spectrum.csv"

## Import the data:

# for example
X_list <- lapply(fileName_list, read.table)
lambda = 1
iterMax = 50

source("PFA_algorithms.R")
Y <- algorithm_4(X_list, lambda, iterMax)










