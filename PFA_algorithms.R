#All the algorithm of PFA
# article link:


####Algorithm 1####

# INPUT:
#
# d, the number of principal components
# X, one biological data type, see article

algorithm_1 <- function(Xdata, d) {
  #NOTE : the code is voluntary very explicite
  
  X <- as.matrix(Xdata)
  n <- dim(X)[2] #sample number
  ones <- matrix(1, n, n) 
  I <- diag(1, n, n)
  Xbar <- X %*% (I - (1/n)*ones)
  eigs <- eigen(Xbar%*%t(Xbar), TRUE)
  U <- eigs$vectors
  values <- eigs$values #values already ordered
  res <- list(U[1:d], values[1:d])
  return(res)
}

X <- as.matrix(data_list[[3]])
d <- dim(X)[2]
res <- algorithm_1(X, d)


n=d
ones <- matrix(1, n, n) 
I <- diag(1, n, n)
temp <- I - (1/n)%*%ones
class(temp)
Xbar <- X %*% temp


####Algorithm 2####