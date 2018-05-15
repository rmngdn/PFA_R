# PFA_algorithms.R
# author: Romain GUEDON

#### Description ####

# Here is the R version of the Pattern Fusion Analysis (PFA)
# article link: https://doi.org/10.1093/bioinformatics/btx176 

# The variables's names are those of the article 
# most of the time

# We need the library MASS for the pseudo-inverse 
library(MASS)


#### Algorithm 1 ####



algorithm_1 <- function(X_i_data,n) {
  # X_i_data is a table with features as rows and patients/samples as columns 
  # n is the number of samples/patients
  
  # Compute X_i_bar:
  
  X_i <- as.matrix(X_i_data)
  h_i <- dim(X_i)[1] #features number
  X_i_bar <- X_i %*% (diag(1, n, n) - (1/n)*matrix(1, n, n))
  
  # Compute the eigen pairs of X_i_bar %*% t(X_i_bar):
  
  if (h_i >= n) { #if the matrix has too much features
    
    #we compute t(X_i_bar) %*% X_i_bar instead of X_i_bar %*% t(X_i_bar) to drastically reduce the time of calcul
    eigs <- eigen(t(X_i_bar)%*%X_i_bar, TRUE)
    
    #init a new data frame with the correct size:
    tmp_vectors <- data.frame(matrix(0, h_i, n))
    
    # the non-zero eigen values of X_i_bar %*% t(X_i_bar) are the same than t(X_i_bar) %*% X_i_bar
    # and for each non-zero eigen value lamda, the eigen vectors of X_i_bar %*% t(X_i_bar) are X_i_bar %*% eigVect, 
    # with eigVect an eigen vector of t(X_i_bar) %*% X_i_bar associated to lambda 
    
    for (i in 1:n) {
      # for each eigen value, we transform the old eigen vector into the good one as previoulsy described:
      tmp_vectors[, i] <- X_i_bar %*% as.matrix(eigs$vectors[, i])
    }
    # replace:
    eigs$vectors <- tmp_vectors
  } else {
    eigs <- eigen(X_i_bar%*%t(X_i_bar), TRUE)
  }
  
  # Find the min. number of eigen values/vectors d_i:
  
  d_i <- 1
  continue <- TRUE
  sumEigen <- sum(eigs$values)

  while(continue) {#it stops for sure when d_i = n or h (if h < n)
    #eigen values are ordered
    if (sum(eigs$values[1:d_i])/sumEigen >= 0.8) {
      continue <- FALSE
    } else {
      d_i <- d_i + 1
    }
  }
  # Compute Y_i:
  Y_i <-  t(as.matrix(eigs$vectors[,1:d_i])) %*% X_i_bar # Ui_di' * X_i_bar
  
  return(list("Y_i" = Y_i, "d_i" = d_i))
}




#### Algorithm 2 ####

# Algorithm to calculate optimal weight matrix W






algorithm_2 <- function(delta,lambda,M) {
  # we're looking for the maximum m s.t. theta - phi_m > 0
  # so we start with m=M, decreasingly and stop when the condition 
  # is satisfied
  for (m in M:1) {
    theta <- (2*lambda + sum(delta[1:m]))/m
    if (theta - delta[m] > 0) {
      N <- m
      break 
    }
  }
  W <- matrix(0,M,1)
  for (m in 1:N) {
    W[m] <- ((2*lambda + sum(delta[1:m]))/m - delta[m])/(2*lambda) 
  }
  return(W)
}


#### Algorithm 3 ####

# it's the iterative updating process for PFA:



algorithm_3 <- function(Ys_list, ds_list, k, n, lambda, maxIter) {
  # Ys_ds_list is the result of the function algorithm_1
  # k is the number of data types
  # n is the number of patients/samples
  # lambda is the tuning parameter
  
  # Calculate d
  d = min(ds_list)
  
  # Initialize W
  M <- k*n
  W <- (1/M)*matrix(1, M, 1) #same weight for every instance (one instance per data type) of every sample/patient
  
  OldErrorValue <- Inf #init
  
  for (iter in 1:maxIter) {
    # Optimize Y according to formula (10) in main text:
    
      # Calculate Phi
    
    
    Phi <- matrix(0, n, n)
    W_List <- vector('list', k)
    V_List <- vector('list', k)
    for (i in 1:k) {
      W_i <- diag(sqrt(W[((i-1)*n + 1):(i*n)]))
      W_List[[i]] <- W_i #used later to calculate delta
      tmp <- (diag(1,n,n) - (1/n)*matrix(1,n,n)) %*%
              (diag(1,n,n) - ginv(Ys_list[[i]]%*%W_i)%*%(Ys_list[[i]]%*%W_i))  
      V_List[[i]] <- norm(tmp,type = "F")#used later to calculate delta
      Phi <- Phi + (1/V_List[[i]])*(tmp%*%t(tmp))
    }
    
    rm(tmp, W_i)
    
      # Compute the eigens of Phi
    
    eigsPhi <- eigen(x = Phi, symmetric = TRUE)
    
    #Y is composed of the d eigen vectors of the 2nd to (d+1)th smallest 
    #eigen values:
    Y <- eigsPhi$vectors[,(n-1):(n-d)]
    Y <- t(Y) # It's not written in the article but it's done in the matlab script and 
              # we can't do the matrix products with Y if we don't transpose!
    
    # Calculate the value of tr(Y*Phi*Y') + lambda*norm2(W)^2 
    # and see if it has decreased compared to the previous iteration
    NewErrorValue <- sum( diag( Y %*% Phi %*% t(Y) ) ) + lambda * t(W) %*% W 
    
    
    if (NewErrorValue - OldErrorValue < 0) {
        #value is decreasing, so we continue
        OldErrorValue <- NewErrorValue
    } else {
      break
    }
    
    
    #Warn the user if it has not converged :
    
    if (iter==maxIter){print("maxIter reached without convergence")}
    
    
    # Optimizing W according to Algorithm 2
    
    
      # Calculate Delta:
    delta <- vector('numeric', M)
    for (i in 1:k) {
      for (j in 1:n) {
        L_i <- 
        tmp2 <- Y[,j]
               - Y%*%matrix(1,n,1)/n 
               - Y%*%(diag(1,n,n) - (1/n)*matrix(1,n,n))%*%ginv(Ys_list[[i]]%*%W_List[[i]])%*%Ys_list[[i]][,j]
                
        delta[i*j] <-(t(tmp2)%*%tmp2)/V_List[[i]]
      }
    }
    rm(tmp2)
    
    #sort it:
    delta <- delta[order(delta)]
    
    W <- algorithm_2(delta, lambda, M)
   
  }
  
  print("Number of iterations = ")
  print(iter)
  return(Y)
}

    
#### Algorithm 4 ####



algorithm_4 <- function(X_list, lambda, maxIter) {
  # X_list = [X_1, X_2, ..., X_k] i.e. the data tables/matrix
  # lamda is the tuning parameter
  # maxIter is the maximum number of iterations for algorithm 3
  
  k = length(X_list) # number of data types
  n = dim(X_list[[1]])[2] # number of samples
  
  # 1. Computing the local sample-spectrum 
  # of each data type according to Algorithm 1:
  Ys_list <- vector('list', k)
  ds_list <- vector('integer', k)
  for (i in 1:k) {
    tmp <- algorithm_1(X_list[[i]], n)
    Ys_list[[i]] <- tmp[[1]]
    ds_list[[i]] <- tmp[[2]]
  }
  rm(tmp)# remove tmp from workspace
  
  # 2. Optimizing Y according to Algorithm 3:
  Y <- algorithm_3(Ys_list, ds_list, k, n, lambda, maxIter)
  colnames(Y) <- colnames(X_list[[1]]) # samples names
  return(Y)
}



