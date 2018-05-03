#All the algorithm of PFA
# article link:


####Algorithm 1####

# INPUT:
#
# d, the number of principal components
# X, one biological data type, see article

algorithm_1 <- function(X_i_data) {
  
  # Compute X_i_bar:
  
  X_i <- as.matrix(X_i_data)
  n <- dim(X_i)[2] #samples number
  h <- dim(X_i)[1] #features number
  ones <- matrix(1, n, n) 
  I <- diag(1, n, n)
  X_i_bar <- X_i %*% (I - (1/n)*ones)
  
  # Compute the eigen pairs:
  
  if (h >= n) { #if the matrix has too much features
    
    #we compute t(X_i_bar) %*% X_i_bar instead of X_i_bar %*% t(X_i_bar) to drastically reduce the size
    eigs <- eigen(t(X_i_bar)%*%X_i_bar, TRUE)
    
    #init a new data frame with the correct size:
    tmp <- data.frame(matrix(0, h, n))
    
    # the non-zero eigen values of X_i_bar %*% t(X_i_bar) are the same than t(X_i_bar) %*% X_i_bar
    # and for each non-zero eigen value lamda, the eigen vectors of X_i_bar %*% t(X_i_bar) are X_i_bar %*% eigVect, 
    # with eigVect an eigen vector of t(X_i_bar) %*% X_i_bar associated to lambda 
    
    for (i in 1:n) {
      # for each eigen value, we transform the old eigen vector into the good one as previoulsy described:
      tmp[, i] <- X_i_bar %*% as.matrix(eigs$vectors[, i])
    }
    # replace:
    eigs$vectors <- tmp
  } else {
    eigs <- eigen(X_i_bar%*%t(X_i_bar), TRUE)
  }
  
  # Find the min. number of eigen values/vectors d_i:
  
  d_i <- 1
  continue <- TRUE
  sumEigen <- sum(eigs$values)

  while(continue) {
    #eigen values are ordered
    if (sum(eigs$values[1:d_i])/sumEigen >= 0.8) {
      continue <- FALSE
    } else {d_i <- d_i + 1}
  }
  # Compute Y_i:
  Y_i <-  t(as.matrix(eigs$vectors[,1:d_i])) %*% X_i_bar
  colnames(Y_i) <- colnames(X_i)
  return(Y_i)
}




####Algorithm 2####
algorithm_2 <- function(Delta) {
  
  
  
}