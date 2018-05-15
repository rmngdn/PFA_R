# Test of PFA
# article of reference:
# https://doi.org/10.1093/bioinformatics/btx176 

pathToYourWD <- "~/Documents/Local_Work/PFA to R"
setwd(pathToYourWD)
source("PFA_algorithm_v2.R")



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

source("PFA_algorithm_v2.R")
Y <- algorithm_4(X_list, lambda, iterMax)


result <- kmeans(t(Y), 3)

Y_matlab <- read.csv("global_sample_spectrum.csv", header = FALSE)

result2 <- kmeans(t(Y_matlab), 3)
myCluster <- result$cluster
matlabCluster <- result2$cluster

