wrapper_HMM_stat_undir_Dens<-function(adjmat, K, thresh, iter.max, theta_init){
  ## This is the combined final updated and working code for EM undirected case for all K
  
  #################################################
  ## Defining a function to update variational parameters gamma using quadratic program solver. Input:
  ## gamma.curr is a N*K matrix for current gamma estimates, pi.curr is a K*1 vector for current pi
  ## estimates, theta.curr is a T_grid*K array for current theta estimates. Given Network is assumed to be an N*N*T_data array
  
  solve_QP_wrapper <- function(quad_lin_coeff,constraint_matrix,constraint_vector,node_ID,gamma.curr){
    try_QP <- try(gamma_next_vec<-solve.QP((-2*diag(as.vector(quad_lin_coeff[node_ID,,1]))),as.vector(quad_lin_coeff[node_ID,,2]),t(constraint_matrix),constraint_vector)$solution,silent = TRUE)
    if('numeric' %in% class(try_QP)){
      gamma_next_vec <- try_QP
    }
    else{
      gamma_next_vec <- gamma.curr[node_ID,]
      # print("red flag")
      # print(paste0("Node ID for red flag is ", node_ID))
    }
    return(gamma_next_vec)
  }
  
  gamma.update.wrapper<-function(gamma.curr,pi.curr,theta.curr,network,N,K){
    gamma.next<-matrix(NA_real_,N,K)
    constraint_matrix<-matrix(NA_real_,2+K,K)
    constraint_matrix[1,]<-rep(1,K)
    constraint_matrix[2,]<-rep(-1,K)
    constraint_matrix[3:(K+2),]<-diag(K)
    constraint_vector<-c(1,-1,rep(0,K))
    
    quad_lin_coeff <- gamma_update_HMM_stat_undir(gamma=gamma.curr, pi=pi.curr, theta=theta.curr, network=network, N=N, K=K)
    
    for (i in 1:N){
      if(sum(is.nan(quad_lin_coeff[i,,1])) > 0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,1])),1] <- -Inf
      }
      if(sum(is.nan(quad_lin_coeff[i,,2])) > 0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,2])),2] <- Inf
      }
      if(sum(quad_lin_coeff[i,,1] <= -10^8) == 0){
        if(sum(quad_lin_coeff[i,,2] >= 10^8) == 0){
          gamma.next[i,] <- solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i, gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2] >= 10^8) > 0) && (sum(quad_lin_coeff[i,,2] >= 10^8) < K)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2] >= 10^8),2] <- 10^10
          gamma.next[i,] <- solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector, node_ID = i, gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2] >= 10^8)==K){
          gamma.next[i,] <- gamma.curr[i,]
        }
      }else if((sum(quad_lin_coeff[i,,1] <= -10^8) > 0) && (sum(quad_lin_coeff[i,,1] <= -10^8) < K)){
        quad_lin_coeff[i,which(quad_lin_coeff[i,,1] <= -10^8),1] <- -(10^10)
        if(sum(quad_lin_coeff[i,,2] >= 10^8) == 0){
          gamma.next[i,] <- solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2] >= 10^8) > 0) && (sum(quad_lin_coeff[i,,2] >= 10^8) < K)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2] >= 10^8),2] <- 10^10
          gamma.next[i,] <- solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2] >= 10^8) == K){
          gamma.next[i,]<-gamma.curr[i,]
        }
      }else if(sum(quad_lin_coeff[i,,1] <= -10^8) == K){
        gamma.next[i,] <- gamma.curr[i,]
      }
      #print(i)
    }
    # for (i in 1:N){
    #   gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
    # }
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma.next <- t(apply(X = gamma.next,MARGIN = 1, FUN = function(x){
      x_norm <- x/sum(x)
      return(x_norm)
    }))
    return(gamma.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of pi
  pi.update<-function(gamma.curr,N,K){
    pi.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## normalization of pi
    pi.next<-pi.next/sum(pi.next)
    return(pi.next)
  }
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update<-function(theta.curr,pi,gamma,network,N,K){
    theta.next<-rep(NA_real_,K)
    gradient<-grad_HMM_stat_undir(theta=as.vector(theta.curr), gamma=gamma, network=network, N=N, K=K)
    hess<-hess_HMM_stat_undir(theta=as.vector(theta.curr), gamma=gamma, N=N, K=K)
    theta.next<-as.vector(theta.curr)-as.vector(solve(hess)%*%gradient)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,network,K,n_iter,thresh,iter.max){
    N<-dim(network)[1]
    
    ## initializing the arrays for parameters
    gamma<-array(NA_real_,dim=c(N,K,n_iter))
    pi<-matrix(NA_real_,K,n_iter)
    theta<-matrix(NA_real_,K,n_iter)
    gamma[,,1]<-start[[1]]
    pi[,1]<-start[[2]]
    theta[,1]<-start[[3]]
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thresh)&(iter_index<iter.max)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the N*K gamma matrix i.e. variational variational parameters
      gamma[,,iter_index]<-gamma.update.wrapper(gamma.curr=gamma[,,iter_index-1],pi.curr=pi[,iter_index-1], theta.curr=theta[,iter_index-1],network=network,N=N,K=K)
      ## Updating the pi vector
      pi[,iter_index]<-pi.update(gamma.curr=gamma[,,iter_index], N=N, K=K)
      ## Updating the theta matrix
      theta[,iter_index]<-theta.update(theta.curr=theta[,iter_index-1], pi=pi[,iter_index],gamma=gamma[,,iter_index], network=network, N=N, K=K)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_HMM_stat_undir(gamma=gamma[,,iter_index], pi=pi[,iter_index], theta = theta[,iter_index], network=network, N=N, K=K)
      error<-abs(((ELBO_grid.prev-ELBO_grid.curr)/ELBO_grid.prev))
      # print(error)
      # print(iter_index)
      # print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(list('variational parameters'=gamma,'mixture proportions'=pi,'network canonical parameters'=theta))
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining the functions for K=1
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  ## Defining the functions for K=1
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update_K1<-function(theta.curr,network,N){
    gradient<-grad_HMM_stat_undir_K1(theta=theta.curr, network=network, N=N)
    hess<-hess_HMM_stat_undir_K1(theta=theta.curr, N=N)
    theta.next<-theta.curr-(gradient/hess)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator_K1<-function(start,network,K,n_iter,thresh,iter.max){
    N<-dim(network)[1]
    ## initializing the arrays for parameters
    theta<-rep(NA_real_,n_iter)
    theta[1]<-start
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thresh)&(iter_index<iter.max)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the theta matrix
      theta[iter_index]<-theta.update_K1(theta.curr=theta[iter_index-1], network=network, N=N)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_HMM_stat_undir_K1(theta = theta[iter_index], network=network, N=N)
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      # print(iter_index)
      # print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(theta)
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining Model Selection functions based on converged paramters
  
  #################################################
  ## Defining a function to calculate the estimate of complete log likelihood
  comp_loglik<-function(gamma=NA,alpha=NA,theta,network,N,K,cluster_ids_est=NA){
    if(K!=1){
      ## 1st term
      t1<-0
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          cluster_id_i<-cluster_ids_est[i]
          cluster_id_j<-cluster_ids_est[j]
          exp_val<-exp(theta[cluster_id_i]+theta[cluster_id_j])
          t1<-t1+((network[i,j]*(theta[cluster_id_i]+theta[cluster_id_j]))-(log(1+exp_val)))
        }
      }
      ## 2nd term
      t2<-0
      for (i in 1:N){
        t2<-t2+log(alpha[cluster_ids_est[i]])
      }
      comp_val<-t1+t2
    }else if(K==1){
      comp_val<-0
      exp_val<-exp(2*theta)
      log_exp_val<-log(1+exp_val)
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          comp_val<-comp_val+((network[i,j]*(2*theta))-log_exp_val)
        }
      }
    }
    return(comp_val)
  }
  
  #################################################
  ## Defining a function to calculate the integrated classification likelihood
  ICL<-function(gamma=NA,alpha=NA,theta,network,N,K,cluster_ids_est=NA){
    if(K!=1){
      t1<-comp_loglik(gamma = gamma, alpha = alpha, theta = theta,network = network,N = N,K = K, cluster_ids_est = cluster_ids_est)
      t2<-K*log((N*(N-1))/2)
      ICL_val<-t1-t2
    }else if(K==1){
      t1<-comp_loglik(theta = theta,network = network,N = N,K = K)
      t2<-log((N*(N-1))/2)
      ICL_val<-t1-t2
    }
    return(ICL_val)
  }
  
  
  ########################################################################################################
  ## Defining the Rand Index and RASE functions for evaluating the performance of clustering and estimation    respectively
  ## Rand Index function
  RI<-function(cluster_ids_est, cluster_ids_true){
    n=length(cluster_ids_est)
    RI_val=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        RI_val=RI_val+as.numeric((cluster_ids_est[i]==cluster_ids_est[j])==(cluster_ids_true[i]==cluster_ids_true[j]))
      }
    }
    RI_mean=RI_val/(n*(n-1)/2)
    return(RI_mean)
  }
  
  #################################################
  ## RASE functions
  RASE_theta<-function(theta_est, theta_true){
    if(is.vector(theta_est)!=1){
      RASE_val<-sqrt(sum((theta_est-theta_true)^2)/dim(theta_est)[1])
    }else{RASE_val<-sqrt(sum((theta_est-theta_true)^2)/length(theta_est))}
    return(RASE_val)
  }
  
  #########################################################################################################
  ## Defining the parameters
  N<-dim(adjmat)[1] ## Number of nodes from the network
  ## Total time points in the network for which we have data in the form of adjacency matrix
  K<-K ## Defining the number of clusters
  ## Setting initial values using package lda which includes mixed membership stochastic block model (MMSB). Using first time point network to run MMSB. Using the mixed membership result of MMSB as our intial gamma. Next we use inital gamma to find initial pi (mixing proportion). Lastly, for network parameter theta we start with 0 matrix.
  # set.seed((2))
  MMSB_result <- mmsb.collapsed.gibbs.sampler(beta.prior = list(1,1),K = K,network = adjmat,alpha = 1/2,num.iterations = 100)
  gamma.start <- with(MMSB_result, t(document_sums)/colSums(document_sums))
  #gamma.start <- matrix(rep(1/K, N*K), N, K)
  ## There is some chance that some component of initial gamma will be exactly 0 which can cause problem in calculating log(gamma). Therfore adding very small amount (10^(-3)) to exact 0 value and rescaling to have sum to 1.
  for(i in 1:N){
    gamma.start[i, which(gamma.start[i,] == 0)] <- rep(10^(-3), length(which(gamma.start[i,] == 0)))
    gamma.start[i,] <- gamma.start[i,]/sum(gamma.start[i,])
  }
  
  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  start<-list()
  
  if(K==1){
    start<-theta_init
    param<-iterator_K1(start=start, network=adjmat, K=K, n_iter=1000, thresh=thresh, iter.max = iter.max)
  }else{
    start[[1]]<-gamma.start
    start[[2]]<-rep(1/K,K) ## alpha (initial distribution)
    start[[3]]<-theta_init ## theta
    #debug(iterator)
    param<-iterator(start=start, network=adjmat, K=K, n_iter=1000, thresh=thresh, iter.max = iter.max)
  }
  
  #################################################
  ## extracting the coverged parameter values and calculating BIC
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    if(K==1){
      temp<-is.na(param[n_iter])
    }else{temp<-is.na(param[[1]][1,1,n_iter])}
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  
  if(K==1){
    param_converge<-param[n_last]
    names(param_converge) <- c('network canonical paramters')
    ICL_val <- ICL(theta = param_converge, network = adjmat,N = N,K = K)
    # if(sim_indicator==1){
    #   if(K==K_true){
    #     RASE_theta<-RASE_theta(theta_est = param_converge,theta_true = theta_true)
    #     output_list<-list("parameters_converged" = param_converge, "ICL_val" = ICL_val, "RASE_theta_val" = RASE_theta)
    #   }else{
    #     output_list<-list("parameters_converged" = param_converge, "ICL_val" = ICL_val)
    #   }
    # }else{output_list<-list("parameters_converged" = param_converge, "ICL_val" = ICL_val)}
    output_list <- list('coefficients' = param_converge, 'ICL' = ICL_val)
  }else{
    param_converge[[1]]<-param[[1]][,,n_last]
    param_converge[[2]]<-param[[2]][,n_last]
    param_converge[[3]]<-param[[3]][,n_last]
    names(param_converge) <- c('variational parameters', 'mixture proportions', 'network canonical paramters')
    cluster_ids_est<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
      cluster_id<-which.max(param_converge[[1]][x,])
      return(cluster_id)
    }))
    ICL_val<-ICL(gamma = param_converge[[1]], alpha=param_converge[[2]] ,theta = param_converge[[3]],network = adjmat,N = N,K = K, cluster_ids_est = cluster_ids_est)
    # if(sim_indicator==1){
    #   RI_val<-RI(cluster_ids_est = cluster_ids_est,cluster_ids_true = cluster_ids_true)
    #   if(K==K_true){
    #     K_permute_mat<-do.call(rbind,permn(1:K))
    #     RASE_theta_vec<-rep(NA_real_,nrow(K_permute_mat))
    #     RASE_pi_vec<-rep(NA_real_,nrow(K_permute_mat))
    #     for (k in 1:nrow(K_permute_mat)){
    #       theta_est<-param_converge[[3]][K_permute_mat[k,]]
    #       RASE_theta_vec[k]<-RASE_theta(theta_est = theta_est,theta_true = theta_true)
    #     }
    #     permute_true_id_theta<-which.min(RASE_theta_vec)
    #     RASE_theta_val<-RASE_theta_vec[permute_true_id_theta]
    #     output_list <- list("parameters_converged" = param_converge, "Estimated_cluster_IDs" = cluster_ids_est,"ICL_val" = ICL_val,"RI_val" = RI_val, "RASE_theta_val" = RASE_theta_val)
    #   }else{
    #     output_list <- list("parameters_converged" = param_converge, "Estimated_cluster_IDs" = cluster_ids_est,"ICL_val" = ICL_val, "RI_val" = RI_val)
    #   }
    # }else{output_list<-list("parameters_converged" = param_converge, "Estimated_cluster_IDs" = cluster_ids_est,"ICL_val" = ICL_val)}
    output_list<-list('coefficients' = param_converge[[3]], 
                      'probability' = param_converge[[1]], 
                      'clust.labels' = cluster_ids_est,
                      'ICL' = ICL_val)
  }
  return(output_list)
}

# directed
wrapper_HMM_stat_dir_Dens<-function(adjmat, K, thresh, iter.max, theta_init){
  ## This is the combined final updated and working code for EM directed case for all K
  
  solve_QP_wrapper <- function(quad_lin_coeff,constraint_matrix,constraint_vector,node_ID,gamma.curr){
    try_QP <- try(gamma_next_vec<-solve.QP((-2*diag(as.vector(quad_lin_coeff[node_ID,,1]))),as.vector(quad_lin_coeff[node_ID,,2]),t(constraint_matrix),constraint_vector)$solution,silent = TRUE)
    if('numeric' %in% class(try_QP)){
      gamma_next_vec <- try_QP
    }
    else{
      gamma_next_vec <- gamma.curr[node_ID,]
      # print("red flag")
      # print(paste0("Node ID for red flag is ", node_ID))
    }
    return(gamma_next_vec)
  }
  
  gamma.update.wrapper<-function(gamma.curr,pi.curr,theta.curr,network,N,K){
    gamma.next<-matrix(NA_real_,N,K)
    constraint_matrix<-matrix(NA_real_,2+K,K)
    constraint_matrix[1,]<-rep(1,K)
    constraint_matrix[2,]<-rep(-1,K)
    constraint_matrix[3:(K+2),]<-diag(K)
    constraint_vector<-c(1,-1,rep(0,K))
    
    quad_lin_coeff <- gamma_update_HMM_stat_dir(gamma=gamma.curr, pi = pi.curr, theta = theta.curr, network = network, N = N, K = K)
    
    for (i in 1:N){
      if(sum(is.nan(quad_lin_coeff[i,,1])) > 0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,1])),1] <- -Inf
      }
      if(sum(is.nan(quad_lin_coeff[i,,2])) > 0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,2])),2] <- Inf
      }
      if(sum(quad_lin_coeff[i,,1] <= -10^8) == 0){
        if(sum(quad_lin_coeff[i,,2] >= 10^8) == 0){
          gamma.next[i,] <- solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i, gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2] >= 10^8) > 0) && (sum(quad_lin_coeff[i,,2] >= 10^8) < K)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2] >= 10^8),2] <- 10^10
          gamma.next[i,] <- solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector, node_ID = i, gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2] >= 10^8)==K){
          gamma.next[i,] <- gamma.curr[i,]
        }
      }else if((sum(quad_lin_coeff[i,,1] <= -10^8) > 0) && (sum(quad_lin_coeff[i,,1] <= -10^8) < K)){
        quad_lin_coeff[i,which(quad_lin_coeff[i,,1] <= -10^8),1] <- -(10^10)
        if(sum(quad_lin_coeff[i,,2] >= 10^8) == 0){
          gamma.next[i,] <- solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2] >= 10^8) > 0) && (sum(quad_lin_coeff[i,,2] >= 10^8) < K)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2] >= 10^8),2] <- 10^10
          gamma.next[i,] <- solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2] >= 10^8) == K){
          gamma.next[i,]<-gamma.curr[i,]
        }
      }else if(sum(quad_lin_coeff[i,,1] <= -10^8) == K){
        gamma.next[i,] <- gamma.curr[i,]
      }
      #print(i)
    }
    # for (i in 1:N){
    #   gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
    # }
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma.next <- t(apply(X = gamma.next,MARGIN = 1, FUN = function(x){
      x_norm <- x/sum(x)
      return(x_norm)
    }))
    return(gamma.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of pi
  pi.update<-function(gamma.curr,N,K){
    pi.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## normalization of pi
    pi.next<-pi.next/sum(pi.next)
    return(pi.next)
  }
  
  #################################################
  ## Defining the update function of the theta with dimensions T_grid*K*2
  theta.update<-function(theta.curr,gamma,network,N,K){
    gradient_oe<-grad_HMM_stat_dir_oe(theta=theta.curr, gamma=gamma, network=network, N=N, K=K)
    gradient_re<-grad_HMM_stat_dir_re(theta=theta.curr, gamma=gamma, network=network, N=N, K=K)
    gradient<-c(gradient_oe,gradient_re)
    hess_oe<-hess_HMM_stat_dir_oe(theta=theta.curr, gamma=gamma, N=N, K=K)
    hess_re<-hess_HMM_stat_dir_re(theta=theta.curr, gamma=gamma, N=N, K=K)
    hess_oe_re<-hess_HMM_stat_dir_oe_re(theta=theta.curr, gamma=gamma, N=N, K=K)
    hess<-matrix(NA_real_,2*K,2*K)
    hess[1:K,1:K]<-hess_oe
    hess[((K+1):(2*K)),((K+1):(2*K))]<-hess_re
    hess[(1:K),((K+1):(2*K))]<-hess_oe_re
    hess[((K+1):(2*K)),(1:K)]<-t(hess_oe_re)
    theta.next<-matrix(c(theta.curr[,1],theta.curr[,2])-as.vector(solve(hess)%*%gradient),K,2)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,network,K,n_iter,thresh,iter.max){
    N<-dim(network)[1]
    ## initializing the arrays for parameters
    gamma<-array(NA_real_,dim=c(N,K,n_iter))
    pi<-matrix(NA_real_,K,n_iter)
    theta<-array(NA_real_,dim=c(K,2,n_iter))
    gamma[,,1]<-start[[1]]
    pi[,1]<-start[[2]]
    theta[,,1]<-start[[3]]
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thresh)&(iter_index<iter.max)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the N*K gamma matrix i.e. variational variational parameters
      gamma[,,iter_index]<-gamma.update.wrapper(gamma.curr=gamma[,,iter_index-1],pi.curr=pi[,iter_index-1], theta.curr=theta[,,iter_index-1],network=network,N=N,K=K)
      ## Updating the pi vector
      pi[,iter_index]<-pi.update(gamma.curr=gamma[,,iter_index], N=N, K=K)
      ## Updating the theta array
      theta[,,iter_index]<-theta.update(theta.curr=theta[,,iter_index-1],gamma=gamma[,,iter_index], network=network, N=N, K=K)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_HMM_stat_dir(gamma=gamma[,,iter_index], alpha=pi[,iter_index], theta = theta[,,iter_index], network=network, N=N, K=K)
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      # print(error)
      # print(iter_index)
      # print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(list(gamma,pi,theta))
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining the functions for K=1
  
  #################################################
  theta.update_K1<-function(theta.curr,network,N){
    gradient_oe<-grad_HMM_stat_dir_oe_K1(theta=theta.curr, network=network, N=N)
    gradient_re<-grad_HMM_stat_dir_re_K1(theta=theta.curr, network=network, N=N)
    hess_oe<-hess_HMM_stat_dir_oe_K1(theta=theta.curr, N=N)
    hess_re<-hess_HMM_stat_dir_re_K1(theta=theta.curr, N=N)
    hess_oe_re<-hess_HMM_stat_dir_oe_re_K1(theta=theta.curr, N=N)
    gradient<-c(gradient_oe,gradient_re)
    hess_mat<-matrix(c(hess_oe,hess_oe_re,hess_oe_re,hess_re),2,2)
    theta.next<-theta.curr-as.vector(solve(hess_mat)%*%gradient)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator_K1<-function(start,network,K,n_iter,thresh, iter.max){
    N<-dim(network)[1]
    
    ## initializing the arrays for parameters
    theta<-matrix(NA_real_,2,n_iter)
    theta[,1]<-start
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thresh)&(iter_index<iter.max)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the theta matrix
      theta[,iter_index]<-theta.update_K1(theta.curr=theta[,iter_index-1], network=network, N=N)
      # print(theta[,iter_index])
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_HMM_stat_dir_K1(theta = theta[,iter_index], network=network, N=N)
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      # print(iter_index)
      # print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(theta)
  }
  
  #########################################################################################################
  #########################################################################################################
  ## Defining Model Selection functions based on converged paramters
  
  ## Defining a function to calculate the estimate of conditional log likelihood
  comp_loglik<-function(gamma=NA,alpha=NA,pi=NA,theta,network,N,K,cluster_ids_est=NA){
    if(K!=1){
      t1<-0
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          cluster_id_i<-cluster_ids_est[i]
          cluster_id_j<-cluster_ids_est[j]
          exp_val_1=exp(theta[cluster_id_i,1])
          exp_val_2=exp(theta[cluster_id_j,1])
          exp_val_3=exp(theta[cluster_id_i,2]+theta[cluster_id_j,2])
          indicator_10=(network[i,j]==1)&(network[j,i]==0)
          indicator_01=(network[i,j]==0)&(network[j,i]==1)
          indicator_11=(network[i,j]==1)&(network[j,i]==1)
          t1<-t1+((indicator_10*theta[cluster_id_i,1])+(indicator_01*theta[cluster_id_j,1])+(indicator_11*(theta[cluster_id_i,2]+theta[cluster_id_j,2]))-log(1+exp_val_1+exp_val_2+exp_val_3))
        }
      }
      t2<-0
      for (i in 1:N){
        t2<-t2+log(alpha[cluster_ids_est[i]])
      }
      comp_val<-t1+t2
    }else if(K==1){
      comp_val<-0
      exp_val_1=exp(theta[1])
      exp_val_2=exp(2*theta[2])
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          indicator_10=(network[i,j]==1)&(network[j,i]==0)
          indicator_01=(network[i,j]==0)&(network[j,i]==1)
          indicator_11=(network[i,j]==1)&(network[j,i]==1)
          comp_val<-comp_val+((indicator_10*theta[1])+(indicator_01*theta[1])+(indicator_11*(2*theta[2]))-log(1+2*exp_val_1+exp_val_2))
        }
      }
    }
    return(comp_val)
  }
  
  ## Defining a function to calculate the integrated classification likelihood
  ICL<-function(gamma=NA,alpha=NA,pi=NA,theta,network,N,K,cluster_ids_est=NA){
    if(K!=1){
      t1<-comp_loglik(gamma = gamma, alpha = alpha, pi = pi, theta = theta,network = network,N = N,K = K,  cluster_ids_est = cluster_ids_est)
      t2<-2*K*log((N*(N-1))/2)
      ICL_val<-t1-t2
    }else if(K==1){
      t1<-comp_loglik(gamma = gamma, alpha = alpha, pi = pi, theta = theta,network = network,N = N,K = K, cluster_ids_est = cluster_ids_est)
      t2<-2*log((N*(N-1))/2)
      ICL_val<-t1-t2
    }
    return(ICL_val)
  }
  
  ########################################################################################################
  ## Defining the Rand Index and RASE functions for evaluating the performance of clustering and estimation    respectively
  ## Rand Index function
  RI<-function(cluster_ids_est, cluster_ids_true){
    n=length(cluster_ids_est)
    RI_val=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        RI_val=RI_val+as.numeric((cluster_ids_est[i]==cluster_ids_est[j])==(cluster_ids_true[i]==cluster_ids_true[j]))
      }
    }
    RI_mean=RI_val/(n*(n-1)/2)
    return(RI_mean)
  }
  
  #################################################
  ## RASE functions
  RASE_theta<-function(theta_est, theta_true){
    if(is.vector(theta_est)!=1){
      RASE_val<-sqrt(sum((theta_est-theta_true)^2)/dim(theta_est)[1])
    }else{RASE_val<-sqrt(sum((theta_est-theta_true)^2)/length(theta_est))}
    return(RASE_val)
  }
  
  #########################################################################################################
  ## Defining the parameters
  N<-dim(adjmat)[1] ## Number of nodes from the network
  K<-K ## Defining the number of clusters
  
  ####################################################
  ## Initializing gamma using MMSB
  if(K==1){
    gamma.start <- rep(1,N)
  }else{
    ## Setting initial values using package lda which includes mixed membership stochastic block model (MMSB). Using first time point network to run MMSB. Using the mixed membership result of MMSB as our intial gamma. Next we use inital gamma to find initial pi (mixing proportion). Lastly, for network parameter theta we start with 0 matrix.
    # set.seed((2))
    MMSB_result <- mmsb.collapsed.gibbs.sampler(beta.prior = list(1,1),K = K,network = adjmat,alpha = 1/2,num.iterations = 100)
    gamma.start <- with(MMSB_result, t(document_sums)/colSums(document_sums))
    #gamma.start <- matrix(rep(1/K, N*K), N, K)
    ## There is some chance that some component of initial gamma will be exactly 0 which can cause problem in calculating log(gamma). Therfore adding very small amount (10^(-3)) to exact 0 value and rescaling to have sum to 1. 
    # for(i in 1:N){
    #   gamma.start[i, which(gamma.start[i,] == 0)] <- rep(10^(-3), length(which(gamma.start[i,] == 0)))
    #   gamma.start[i,] <- gamma.start[i,]/sum(gamma.start[i,])
    # }
  }
  
  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  start<-list()
  
  if(K==1){
    start<-theta_init
    param<-iterator_K1(start=start, network=adjmat, K=K, n_iter=1000, thresh=thresh, iter.max = iter.max)
  }else{
    start[[1]]<-gamma.start
    start[[2]]<-rep(1/K,K) ## alpha (initial distribution)
    start[[3]]<-theta_init
    param<-iterator(start=start, network=adjmat, K=K, n_iter=1000, thresh=thresh, iter.max = iter.max)
  }
  
  #################################################
  ## extracting the coverged parameter values and calculating BIC
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    if(K==1){
      temp<-is.na(param[1,n_iter])
    }else{temp<-is.na(param[[1]][1,1,n_iter])}
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  if(K==1){
    param_converge<-param[,n_last]
    ICL_val<-ICL(theta = param_converge,network = adjmat,N = N,K=K)
    # if(sim_indicator==1){
    #   if(K==K_true){
    #     RASE_theta<-RASE_theta(theta_est = param_converge,theta_true = theta_true)
    #     output_list<-list(param_converge,ICL_val,RASE_theta)
    #   }else{
    #     output_list<-list(param_converge,ICL_val)
    #   }
    # }else{output_list<-list(param_converge,ICL_val)}
    output_list<-list('coefficients' = param_converge, 'ICL' = ICL_val)
  }else{
    param_converge[[1]]<-param[[1]][,,n_last]
    param_converge[[2]]<-param[[2]][,n_last]
    param_converge[[3]]<-param[[3]][,,n_last]
    cluster_ids_est<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
      cluster_id<-which.max(param_converge[[1]][x,])
      return(cluster_id)
    }))
    
    ICL_val<-ICL(gamma = param_converge[[1]], alpha=param_converge[[2]] ,theta = param_converge[[3]],network = adjmat,N = N,K = K, cluster_ids_est = cluster_ids_est)
    # if(sim_indicator==1){
    #   RI_val<-RI(cluster_ids_est = cluster_ids_est,cluster_ids_true = cluster_ids_true)
    #   if(K==K_true){
    #     K_permute_mat<-do.call(rbind,permn(1:K))
    #     RASE_theta_vec<-rep(NA_real_,nrow(K_permute_mat))
    #     for (k in 1:nrow(K_permute_mat)){
    #       theta_est<-param_converge[[3]][K_permute_mat[k,],]
    #       RASE_theta_vec[k]<-RASE_theta(theta_est = theta_est,theta_true = theta_true)
    #     }
    #     permute_true_id_theta<-which.min(RASE_theta_vec)
    #     RASE_theta_val<-RASE_theta_vec[permute_true_id_theta]
    #     output_list<-list(param_converge,cluster_ids_est,ICL_val,RI_val,RASE_theta_val)
    #   }else{
    #     output_list<-list(param_converge,cluster_ids_est,ICL_val,RI_val)
    #   }
    # }else{output_list<-list(param_converge,cluster_ids_est,ICL_val)}
    output_list<-list('coefficients' = param_converge[[3]], 
                      'probability' = param_converge[[1]], 
                      'clust.labels' = cluster_ids_est,
                      'ICL' = ICL_val)
  }
  return(output_list)
}

############ Main Functions ############
# ergmclust() 
ergmclust <- function(adjmat,K,directed=FALSE,thresh=1e-06,iter.max=200,coef.init=NULL){
  if (directed){
    if(is.null(coef.init)){
      wrapper_HMM_stat_dir_Dens(adjmat=adjmat,K=K,thresh=thresh, iter.max=iter.max, theta_init = jitter(matrix(0,K,2)))
    }else{
      wrapper_HMM_stat_dir_Dens(adjmat=adjmat,K=K,thresh=thresh, iter.max=iter.max, theta_init = coef.init)
    }
  }
  else {
    if (is.null(coef.init)){
      wrapper_HMM_stat_undir_Dens(adjmat=adjmat,K=K,thresh=thresh,iter.max=iter.max,theta_init = jitter(rep(0,K)))
    }
    else {
      wrapper_HMM_stat_undir_Dens(adjmat=adjmat,K=K,thresh=thresh,iter.max=iter.max,theta_init = coef.init)
    }
  }
}

# graphics tool
ergmclust.plot <- function(adjmat, K, directed=FALSE, thresh=1e-06, iter.max=200, coef.init=NULL, node.labels = NULL){
  net_result <- list()
  network_object <- make_empty_graph()
  n_clust <- numeric()
  if (directed){
    if(is.null(coef.init)) {
      net_result<- ergmclust(adjmat=adjmat, K=K, directed=T, thresh=thresh, iter.max=iter.max)
      network_object <- graph.adjacency(adjmat, mode = 'directed')
    }
    else {
      net_result<- ergmclust(adjmat=adjmat, K=K, directed=T, thresh=thresh, iter.max=iter.max, coef.init=coef.init)
      network_object <- graph.adjacency(adjmat, mode = 'directed')
    }
  }
  else {
    if(is.null(coef.init)) {
      net_result <- ergmclust(adjmat=adjmat, K=K, directed=F, thresh=thresh, iter.max=iter.max, coef.init=NULL)
      network_object <- graph.adjacency(adjmat, mode = 'undirected')
    }
    else {
      net_result <- ergmclust(adjmat=adjmat, K=K, directed=F, thresh=thresh, iter.max=iter.max, coef.init=coef.init)
      network_object <- graph.adjacency(adjmat, mode = 'undirected')
    }
  }
  
  V(network_object)$label.cex <- 0.6
  V(network_object)$label.font <- 2
  V(network_object)$Est = net_result[["clust.labels"]]
  n_clust=length(unique(net_result[["clust.labels"]]))
  
  for (i in 1:dim(adjmat)[1]){
    for (k in 1:n_clust){
      if (V(network_object)$Est[i] == k) {
        V(network_object)$color[i] = rainbow(n_clust)[k]
      }
    }
  }
  if(is.null(node.labels)){
    plot(network_object, vertex.frame.color=NA, 
         vertex.size=7,
         edge.width=0.5, 
         edge.arrow.size=0.15,
         vertex.label=1:dim(adjmat)[1], 
         vertex.label.cex=0.8,
         #vertex.label=NA, 
         layout=layout.fruchterman.reingold(network_object))
  }
  else {
    plot(network_object, vertex.frame.color=NA, 
         vertex.size=7,
         edge.width=0.5, 
         edge.arrow.size=0.05,
         vertex.label=node.labels,
         vertex.label.cex=0.8,
         layout=layout.fruchterman.reingold(network_object))
  }
}

# model selection based on ICL
ergmclust.ICL <- function(adjmat, Kmax=5, directed=FALSE, thresh=1e-06, iter.max=200, coef.init = NULL){
  ICL.vec <- numeric(Kmax) # a vector store the ICL values
  net.res <- list()
  for (K_index in 1:Kmax){
    net.res[[K_index]] <- ergmclust(adjmat=adjmat, K=K_index, directed=directed, thresh=thresh, iter.max=iter.max, coef.init=coef.init)
    ICL.vec[K_index] <- net.res[[K_index]][['ICL']]
  }
  chosen.clust <- which.max(ICL.vec)
  net_result <- list('Kselect'=chosen.clust, 
                     'coefficients'=net.res[[chosen.clust]][['coefficients']],
                     'probability'=net.res[[chosen.clust]][['probability']],
                     'clust.labels'=net.res[[chosen.clust]][['clust.labels']],
                     'ICL'=net.res[[chosen.clust]][['ICL']]
  )
  return(net_result)
}