# PFA_algorithms.R
# author: Romain GUEDON

#### Description ####

# Here is the R version of the Pattern Fusion Analysis (PFA)
# article link:

# The variables's names are those of the article 
# most of the time


#### Algorithm 1 ####

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
  return(list(Y_i, d_i))
}





#### Algorithm 2 ####

# Algorithm to calculate optimal weight matrix W






# phiCalc <- function(Y_i, )


algorithm_2 <- function(Delta,M) {
  for (N in 1:M) {
    
  }
  
  
}


#### Algorithm 3 ####

algorithm_3 <- function(Ys_ds_list, k, n, lambda) {
  # Ys_ds_list is the result of the function algorithm_1
  # k is the number of data types
  # n is the number of patients/samples
  # lambda is the tuning parameter
  
  # Extract Ys and ds 
  Ys_list <- Ys_ds_list[[1]]
  ds_list <- Ys_ds_list[[2]]
  
  # Set d
  d = min(ds_list)
  
  # Initialize W
  M <- k*n
  W <- matrix(1, M, 1)
  
  
  
  # Optimize Y according to formula (10) in main text:
  
    # Calculate Phi
  
  
  Phi <- matrix(0,n,n)
  for (i in 1:k) {
    W_i <- diag(W[(i-1)*n + 1, i*n])
    tmp <- (diag(1,n,n) - (1/n)*matrix(1,n,n)) 
            %*%
            (diag(1,n,n) - ginv(Ys_list[[i]]%*%W_i)%*%(Ys_list[[i]]%*%W_i))  
    
    Phi <- Phi + (1/norm(tmp,type = "F"))*(tmp%*%t(tmp))
  }
  
    # Compute the eigens of Phi
  
  eigsPhi <- eigen(x = Phi, symmetric = TRUE)
  Y <- eigsPhi$vectors[n:n-d+1] #vectors of the d smallest eigen values
  
  
  # Optimizing W according to Algorithm 2
  
  

  
  
}

    
#### Algorithm 4 ####

# it's the terative updating process for PFA

algorithm_4 <- function(X_list, d_list, lambda) {
  # X_list = [X_1, X_2, ..., X_k]
  # d_list = [d_1, d_2, ..., d_k]
  # lamda is the tuning parameter
  
  # Computing the local sample-spectrum 
  # of each data type according to Algorithm 1:
  local_Ys <- lapply(X_list, algorithm_1)
  
  
  # Optimizing Y according to Algorithm 3:
  algorithm_3()
  
  
}