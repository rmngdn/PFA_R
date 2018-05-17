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
Y_PFA <- algorithm_4(X_list, lambda, iterMax)


#set the function:
clustFunc <- kmeans
res <- clusGap(t(Y_PFA), clustFunc, 5, B = 200, d.power = 2, spaceH0 = "scaledPCA")       
plot(res)
print(res)

# Extract the optimal number of cluster

K_PFA <- 3 # see the print for the good value

#Cluster
result_PFA <- clustFunc(t(Y_PFA), K_PFA)
cluster_PFA <- result_PFA$cluster

# Let's plot the patients against the 2 biggest principal components:

forPlot <- data.frame(t(Y_PFA))
ggplot(forPlot) + aes(x=forPlot[,1], y = forPlot[,2]) + geom_point() 
