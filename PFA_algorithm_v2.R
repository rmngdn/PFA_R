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
  
  X_i <- as.matrix(X_i_data)
  h_i <- dim(X_i)[1] #features number
  
  # Centralize X_i for each feature ? 
  
  #X_mean <- matrix(rowMeans(X_i)) %*% matrix(1,1,n)
  #X_i <- X_i - X_mean
  
  #
   
  # Compute the eigen pairs 
  
  if (h_i >= n) { #if the matrix has too much features
    
    R = t(X_i)%*%X_i/(n-1) 
    
    eigs <- eigen(R, TRUE)
    
    #init a new data frame with the correct size:
    tmp_vectors <- matrix(0, h_i, n)
    #View(eigs$values)
    
    # for each eigen value, we transform the old eigen vector into the good one as previoulsy described:
    for (j in 1:n) {
      
  
      if (eigs$values[j]>=0) {
        # we avoid false negative values for the square root: (artefact of computing, this is eigen values close to machine epsilon)
        # this values won't be useful anyway cause we'll take the eigen vectors associated to the biggest eigen values...
        tmp_vectors[, j] <- (1/sqrt((n-1)*eigs$values[j])) * X_i %*% eigs$vectors[, j]
      }
    }
    # replace:
    eigs$vectors <- tmp_vectors
  } else {
    eigs <- eigen(X_i%*%t(X_i), TRUE)
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
  View(eigs$vectors[1:10])
  # Compute Y_i: #HAVE TO REMOVE THE MINUS SIGN after tests 
  Y_i <-  -t(as.matrix(eigs$vectors[,1:d_i])) %*% X_i # Ui_di' * X_i {we put a minus sign to have the same thing than matlab but it doesnt change anything nornally}
  return(list("Y_i" = Y_i, "d_i" = d_i))
}




#### Algorithm 2 ####

# Algorithm to calculate optimal weight matrix W






algorithm_2 <- function(delta, permutation,lambda, M) {
  # we're looking for the maximum m s.t. theta - phi_m > 0
  # so we start with m=M, decreasingly and stop when the condition 
  # is satisfied
  N <- Inf
  for (m in M:1) {
    theta <- (2*lambda + sum(delta[1:m]))/m
    if (theta - delta[m] >= 0) {
      N <- m
      break 
    }
  }
  W <- matrix(0,M,1)
  W_tmp <- matrix(0, M , 1)
  theta <- (2*lambda + sum(delta[1:N]))/N
  for (m in 1:N) {
    W_tmp[m] <- (theta - delta[m])/(2*lambda) 
  }
  # Then we place the values to the corresponding place in W
  W[permutation] <- W_tmp
  View(W)
  return(W)
}


#### Algorithm 3 ####

# it's the iterative updating process for PFA:



algorithm_3 <- function(Y_list, d_list, k, n, lambda, maxIter) {
  # Ys_d_list is the result of the function algorithm_1
  # k is the number of data types
  # n is the number of patients/samples
  # lambda is the tuning parameter
  
  # Calculate d
  d = min(d_list)
  
  # Initialize W
  M <- k*n
  W <- (1/M)*matrix(1, M, 1) #same weight for every instance (one instance per data type) of every sample/patient
  
  final_err <- Inf #init
  
  
  for (iter in 1:maxIter) {
    # Optimize Y according to formula (10) in main text:
    
    # Calculate Phi
    
    #init Phi, list of S_i and list of W_i:
    Phi <- matrix(0, n, n) # M in matlab code 
    W_List <- vector('list', k)
    S_List <- vector('list', k)
    
    for (i in 1:k) {
      # Extract W_i:
      W_i <- diag(sqrt(W[((i-1)*n + 1):(i*n)]), n, n)
      W_List[[i]] <- W_i #used later to calculate delta
      
      # Compute
      
      #Let's decompose it :
      S_List[[i]] <- sum(diag( (diag(1,n,n) - ginv(Y_list[[i]])%*%Y_list[[i]])%*%t(diag(1,n,n) - ginv(Y_list[[i]])%*%Y_list[[i]])))
      print("S_i = ")
      print(S_List[[i]])
      tmp2 <- ginv(Y_list[[i]]%*%W_i) %*% (Y_list[[i]]%*%W_i)
      tmp <- W_i*(diag(1,n,n) - ginv(Y_list[[i]]%*%W_i)%*%(Y_list[[i]]%*%W_i) )
      Phi <- Phi + (1/S_List[[i]])*tmp%*%t(tmp)
      print("tmp2 = ")
      print(tmp2[1:5,1:5])
      rm(tmp2)
    }
    View(W_List[[1]])
    #View(S_List)
    
    #View(Phi)
    rm(tmp, W_i)
    
    # Compute the smallest eigen values (and associated vector) of Phi
    
    eigsPhi <- eigen(x = Phi, symmetric = TRUE)
    #View(eigsPhi$values)
    #View(eigsPhi$vectors)
    #Y is composed of the d eigen vectors of the 2nd to (d+1)th smallest 
    #eigen values:
    Y <- eigsPhi$vectors[,(n-1):(n-d)]
    Y <- t(Y) # It's not written in the article but it's done in the matlab script and 
              # we can't do the matrix products with Y if we don't transpose!
    #View(Y)
    # Calculate the value of tr(Y*Phi*Y') + lambda*norm2(W)^2 
    # and see if it has decreased compared to the previous iteration
    err <- sum( diag( Y %*% Phi %*% t(Y) ) ) 
    
    
    if (err - final_err < 0) {
      #value is decreasing, so we continue
      final_err <- err
    } else {
      break
    }
    
    
  # Optimizing W according to Algorithm 2
    
    
    # Calculate Delta:
    delta <- vector('numeric', M)
    for (i in 1:k) {
      L_i <- ( Y%*%W_List[[i]] ) %*% ginv( Y_list[[i]]%*%W_List[[i]] )
      for (j in 1:n) {
        tmp2 <- Y[,j] - L_i%*%Y_list[[i]][, j]
        # we sort the values like this: jth patient in ith data is in ((i-1)*n + j)th position
        delta[(i-1)*n + j] <- (1/S_List[[i]])*W[(i-1)*n + j]%*%t(tmp2)%*%tmp2
      }
    }
    rm(tmp2)
    
    #sort it:
    permutation <- order(delta) #keep this permutation 
    delta <- delta[permutation] # ascending order
    
    # Update W:
    W <- algorithm_2(delta, permutation,lambda, M)
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
  Y_list <- vector('list', k)
  d_list <- vector('integer', k)
  for (i in 1:k) {
    tmp <- algorithm_1(X_list[[i]], n)
    Y_list[[i]] <- tmp[[1]]
    d_list[[i]] <- tmp[[2]]
  }
  
  rm(tmp)# remove tmp from workspace
  
  # 2. Optimizing Y according to Algorithm 3:
  Y <- algorithm_3(Y_list, d_list, k, n, lambda, maxIter)
  colnames(Y) <- colnames(X_list[[1]]) # samples names
  return(Y)
}



