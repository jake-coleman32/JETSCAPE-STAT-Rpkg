## Metropolis sampler to carry out calibration ##
calibration_sampler <- function(niters = 1E4,
                                burnin_prop = 0.3,
                                t_kap = 0.1,
                                
                                design_output = NULL,
                                GP_list = NULL, #Output from train_GPs
                                y_exp = NULL,
                                Sigma_E = NULL,
                                
                                include_cov_extra = TRUE,
                                proposal_cor_mat = diag(dim(GP_list[[1]]@input)[2]),
                                upper_theta_trunc = 1,
                                bad_corners = list()){
  
  
  #Renaming for ease of use
  Y <- design_output
  
  y_sd <- apply(Y, 2, sd)
  y_means <- apply(Y, 2, mean)
  
  Y_final <- as.matrix(sweep(Y, 2, y_means, FUN = '-')) %*% diag(1/y_sd)
  Y_svd <- svd(Y_final)
  
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  
  R <- length(GP_list)
  
  #Do various sampler setups
  time_start <- proc.time()
  
  burnin <- burnin_prop*niters
  
  t_out <- matrix(0,niters,dim(proposal_cor_mat)[2])
  t_ratio <- ob_ratio <- 0
  
  #draw from priors
  (t_cur <- t(runif(dim(proposal_cor_mat)[2],0,upper_theta_trunc)))
  
  #Get current values
  (X_cur <- t_cur)
  
  
  #Account for extra variation lost. Remember Y is a (n x q) matrix
  if(include_cov_extra){
    if(R < q){
      V_b = Y_svd$v[ ,(R+1):q]
      S_b = diag(Y_svd$d[(R+1):q])
      cov_extra <- 1/(n-1) * diag(y_sd) %*% V_b %*% S_b^2 %*% t(V_b) %*% diag(y_sd)
    }else{
      cov_extra = 0
    }
  }else{
    cov_extra = 0
  }
  
  
  cur_preds = exp_predictions(theta = X_cur,
                              GP_list = GP_list,
                              Y_svd = Y_svd,
                              y_sd = y_sd,
                              y_means = y_means,
                              cov_extra = cov_extra)
  
  for(i in 1:(niters + burnin)){
    #Get the time every ten percent of iterations
    if(!i%%(niters*0.1)){
      flush.console()
      cat("\r i = ", i, "Elapsed Time:",
          as.numeric((proc.time() - time_start)[3]),
          ' Acc Rate:', t_ratio/i
      ) 
    }
    
    #Propose new t_st
    (t_st <- mvtnorm::rmvnorm(1,mean=as.numeric(t_cur),sigma = t_kap*proposal_cor_mat))
    
    #Check to make sure parameter is within [0,1]
    auto_reject = 0
    if(sum(t_st < 0) + sum(t_st > upper_theta_trunc)){
      auto_reject = 1
    }
    
    #If you need to cut off a corner because MCMC gets stuck there
    if(length(bad_corners)){
      for(corner in bad_corners){
        
        #Get the pair for the corner you care about
        (param_index <- which(!is.na(corner[[1]])))
        if(length(param_index)!=2){
          stop('bad corner must be pairwise')
        }
        
        #get direction and location of the bad corner
        (direction = corner[[1]][param_index])
        (loc = corner[[2]][param_index])
        (corner_test <- t_st[param_index])
        if(sum(is.na(direction)) + sum(is.na(loc))){
          stop('NA in bad corner somewhere')
        }
        
        #See if the proposed points lie in the corner
        if(prod(direction==c('low','low'))){
          auto_reject = auto_reject + (corner_test[1] < loc[1] & corner_test[2] < loc[2])
        }else if(prod(direction==c('low','high'))){
          auto_reject = auto_reject + (corner_test[1] < loc[1] & corner_test[2] > loc[2])
        }else if(prod(direction == c('high','low'))){
          auto_reject = auto_reject +  (corner_test[1] > loc[1] & corner_test[2] < loc[2])
        }else if(prod(direction == c('high','high'))){
          auto_reject = auto_reject + (corner_test[1] > loc[1] & corner_test[2] > loc[2])
        }
        
      }
    }
    if(auto_reject){
      ob_ratio = ob_ratio + 1
    }
    
    
    ### Propose a new calibration parameter ###
    X_st <- t_st
    
    ## Get prediction (mean/uncertainty) at new point
    st_preds = exp_predictions(theta = X_st,
                               GP_list = GP_list,
                               Y_svd = Y_svd,
                               y_sd = y_sd,
                               y_means = y_means,
                               cov_extra = cov_extra)
    
    y_st_means <- st_preds$mean
    y_st_cov <- st_preds$cov
    
    y_cur_means <- cur_preds$mean
    y_cur_cov <- cur_preds$cov
    
    
    #Get likelihoods
    (like_st <- mvtnorm::dmvnorm(x = as.numeric(y_exp),
                                 mean = as.numeric(y_st_means),
                                 sigma = Sigma_E + y_st_cov, #cov_extra_phys_cal is part of y_st_cov, and y_cur_cov
                                 log = TRUE))
    (like_cur <- mvtnorm::dmvnorm(x = as.numeric(y_exp),
                                  mean = as.numeric(y_cur_means),
                                  sigma = Sigma_E + y_cur_cov,
                                  log = TRUE))
    
    
    #Includes constant prior 
    (ratio <- sum(like_st) - sum(like_cur))
    
    if(rexp(1) > -ratio & !auto_reject){
      cur_preds <- st_preds
      t_cur <- t_st
      t_ratio <- t_ratio + 1
    }
    
    if(i > burnin){
      t_out[i-burnin,] <- t_cur
    }
  }#i loop
  
  print("t_ratio is")
  print(t_ratio/(niters + burnin))
  print("ob ratio is")
  print(ob_ratio/(niters + burnin))
  
  return(t_out)
  
}
