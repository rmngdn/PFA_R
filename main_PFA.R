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
data_list <- lapply(fileName_list, read.table)
## Capture the local sample-spectrum for each biological data type

#local_Ys <- lapply(data_list, algorithm_1)

Y2 <- algorithm_1(data_list[[2]])
Y3 <- algorithm_1(data_list[[3]])
#Y1 <- algorithm_1(data_list[[1]]) #too big

## WHY IS IT SO LONG ? ##

X <- as.matrix(data_list[[1]])
n <- dim(X)[2]
X_bar <- X %*% (diag(1,n,n) - (1/n)*matrix(1,n,n))
R <- X_bar %*% t(X_bar) #1.1Gb

