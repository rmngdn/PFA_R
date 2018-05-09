# PFA_algorithms.R
# author: Romain GUEDON

#### Description ####

# Here is the R version of the Pattern Fusion Analysis (PFA)
# article link:

# The variables's names are those of the article 
# most of the time


#### Algorithm 1 ####



algorithm_1 <- function(X_i_data) {
  
  # Compute X_i_bar:
  
  X_i <- as.matrix(X_i_data)
  n <- dim(X_i)[2] #samples number
  h <- dim(X_i)[1] #features number
  ones <- matrix(1, n, n) # 11' with 1' = (1,1,1...,1) [1,n]
  I <- diag(1, n, n)
  X_i_bar <- X_i %*% (I - (1/n)*ones)
  
  # Compute the eigen pairs:
  
  if (h >= n) { #if the matrix has too much features
    
    #we compute t(X_i_bar) %*% X_i_bar instead of X_i_bar %*% t(X_i_bar) to drastically reduce the size
    eigs <- eigen(t(X_i_bar)%*%X_i_bar, TRUE)
    
    #init a new data frame with the correct size:
    tmp_vectors <- data.frame(matrix(0, h, n))
    
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

  while(continue) {#it stops for sure when d_i = n
    #eigen values are ordered
    if (sum(eigs$values[1:d_i])/sumEigen >= 0.8) {
      continue <- FALSE
    } else {
      d_i <- d_i + 1
    }
  }
  # Compute Y_i:
  Y_i <-  t(as.matrix(eigs$vectors[,1:d_i])) %*% X_i_bar # Ui_di' * X_i_bar
  colnames(Y_i) <- colnames(X_i)
  return(list(Y_i, d_i))
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



algorithm_3 <- function(Ys_ds_list, k, n, lambda, maxIter) {
  # Ys_ds_list is the result of the function algorithm_1
  # k is the number of data types
  # n is the number of patients/samples
  # lambda is the tuning parameter
  
  # Extract Ys 
  Ys_list <- Ys_ds_list[[1]]
  
  # Calculate d
  d = min(Ys_ds_list[[2]])
  
  # Initialize W
  M <- k*n
  W <- (1/M)*matrix(1, M, 1)
  
  
  for (iter in 1:maxIter) {
    # Optimize Y according to formula (10) in main text:
    
      # Calculate Phi
    
    
    Phi <- matrix(0, n, n)
    W_List <- vector('list', k)
    V_List <- vector('list', k)
    for (i in 1:k) {
      W_i <- diag(W[(i-1)*n + 1, i*n])
      W_List[i] <- W_i #used later to calculate delta
      tmp <- (diag(1,n,n) - (1/n)*matrix(1,n,n)) 
              %*%
              (diag(1,n,n) - ginv(Ys_list[[i]]%*%W_i)%*%(Ys_list[[i]]%*%W_i))  
      V_List[i] <- norm(tmp,type = "F")#used later to calculate delta
      Phi <- Phi + (1/V_List[i])*(tmp%*%t(tmp))
    }
    
      # Compute the eigens of Phi
    
    eigsPhi <- eigen(x = Phi, symmetric = TRUE)
    #Y is composed of the d eigen vectors of the 2nd to (d+1)th smallest 
    #eigen values:
    Y <- eigsPhi$vectors[(n-1):(n-d)] 
    
    
    # Calculate the value of tr(Y*Phi*Y') + lambda*norm2(W)^2 
    # and see if it has decreased compared to the previous iteration
    
    NewObjectiveValue <- sum(diag(Y%*%Phi%*%t(Y))) + lamda*t(W)%*%W 
    
    if (iter==1){
      OldObjectiveValue <- NewObjectiveValue
    } else {
      if (NewObjectiveValue - OldObjectiveValue < 0) {
        #value is decreasing, so we continue
        OldObjectiveValue <- NewObjectiveValue
      } else {
        break 
      }
    }
    
    #Warn the user if it has not converged :
    
    if (iter==maxIter){print("maxIter reached without convergence")}
    
    # Optimizing W according to Algorithm 2
    
    
      # Calculate Delta:
    delta <- vector('list', M)
    for (i in 1:k) {
      for (j in 1:n) {
        delta[i*j] <-(norm(Y[,j]
                          - Y%*%matrix(1,n,1)
                          - Y%*%
                            (diag(1,n,n) - (1/n)*matrix(1,n,n))
                             %*%
                             ginv(Ys_list[[i]]%*%W_List[[i]])
                             %*%
                             Ys_list[[i]][,j]
                          , "F")^2)/norm(V_List[[i]],"F")
      }
    }
    #sort it:
    delta <- delta[order(delta)]
    
    W <- algorithm_2(delta, lambda, M)
  }
  return(Y)
}

    
#### Algorithm 4 ####



algorithm_4 <- function(X_list, lambda, iterMax, k, n) {
  # X_list = [X_1, X_2, ..., X_k] i.e. the data tables/matrix
  # lamda is the tuning parameter
  # iterMax is the maximum number of iterations for algorithm 3
  
  # Computing the local sample-spectrum 
  # of each data type according to Algorithm 1:
  Ys_ds_list <- lapply(X_list, algorithm_1)
  k = length(Ys_ds_list)
  n = Ys_ds_list[[1]]
  
  # Optimizing Y according to Algorithm 3:
  Y <- algorithm_3(Ys_ds_list, k, n, lambda, maxIter)
}