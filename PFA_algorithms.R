#All the algorithm of PFA
# article link:


####Algorithm 1####

# INPUT:
#
# d, the number of principal components
# X, one biological data type, see article

algorithm_1 <- function(X_i_dat) {
  
  # NOTE : the code is voluntary very explicite
  
  X_i <- as.matrix(X_i_dat)
  n <- dim(X_i)[2] #sample number
  ones <- matrix(1, n, n) 
  I <- diag(1, n, n)
  X_i_bar <- X_i %*% (I - (1/n)*ones)
  eigs <- eigen(X_i_bar%*%t(X_i_bar), TRUE)
  
  #let's find the min. number of eigen values/vectors d_i
  d_i <- 1
  continue <- TRUE
  sumEigen <- sum(eigs$values)
  return(eigs)
  # while(d_i <= n && continue) {
  #   #eigen values are ordered
  #   if (sum(eigs$values[1:d_i])/sumEigen >= 0.8) {
  #     continue <- FALSE
  #   } else {d_i <- d_i + 1}
  # }
  # if (!continue) { #condition was satisfied
  #   #compute Y_i
  #   #tmp <- as.matrix(eigs$vectors[,1:d_i])
  #   Y_i <-  t(as.matrix(eigs$vectors[,1:d_i])) %*% X_i_bar
  #   return(Y_i)
  # 
  # } else {
  #   print("Dimensions were not reduced")
  #   # We still return a result ?
  #   Y_i <- t(as.matrix(eigs$vectors)) %*% X_i_bar
  #   return(Y_i)
  # }
}

## TESTS ## 
# X <- as.matrix(data_list[[3]])
# d <- dim(X)[2]
# res <- algorithm_1(X)



####Algorithm 2####