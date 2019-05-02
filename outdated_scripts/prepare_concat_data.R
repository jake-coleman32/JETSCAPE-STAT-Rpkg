library(MASS)
library(cluster)
library(emdbook)
library(dplyr)
library(reshape2)
library(RobustGaSP)
library(mvtnorm)
library(Matrix)
library(stringr)


concatenate_dsets <- function(folder,
                              dsets,
                              design_file,
                              range_file,
                              write_csvs = FALSE,
                              design_names = NULL,
                              est_err_ell = TRUE,
                              ell_err = NULL,
                              subset_high_pT_pbpb = FALSE,
                              errors_are_sd = TRUE,
                              add_header = FALSE){
  
  (design <- read.table(design_file,header = TRUE))
  load(range_file)
  
  scaled_d <- matrix(0,dim(design)[1],dim(design)[2]) %>%
    as.data.frame()
  for(j in 1:dim(design)[2]){
    scaled_d[,j] <- (design[,j] - ranges[[j]][1])/(ranges[[j]][2] - ranges[[j]][1])
  }
  if(!is.null(names(ranges))){
    colnames(scaled_d) <- names(ranges)
  }else if(!is.null(design_names)){
    print(paste0('Using design_names: ',paste0(design_names,collapse = ', ')))
    print('Note: these are also the names for ranges')
    colnames(scaled_d) <- design_names
    names(ranges) <- design_names
  }else{
    print('Warning: No names for scaled_d or ranges')
  }
  
  n_design <- dim(design)[1]
  
  covs <- au <- pT_scaled <- vector('list',length(dsets))
  
  first_concat = TRUE
  for(i in 1:length(dsets)){
    (current_dset = dsets[i])
    (current_experiment = strsplit(current_dset,'-')[[1]][1])
    
    au[[i]] <- read.table(paste0(folder,current_dset,".dat"),header = add_header)
    (names(au[[i]]) <- c("pT","RAA_exp","Stat_err","Sys_err",paste0("RAA_",as.character(1:n_design))))
    
    if(subset_high_pT_pbpb & str_detect(current_experiment,'PbPb')){
      au[[i]] <- au[[i]][which(au[[i]]$pT>30),]
    }
    
    
    if(write_csvs){
      #Write csv of experimental data for easily manipulation later
      write.csv(au[[i]][,1:4],file=paste0(folder,current_dset,"_exp.csv"),row.names = FALSE)
    }
    
    #Separate output, change from wide to long
    mod_dat <- dplyr::select(au[[i]], -c(RAA_exp,Stat_err,Sys_err)) %>%
      melt(id = "pT") %>%
      arrange(by = pT)
    
    (exp_dat <- au[[i]]$RAA_exp)
    
    
    if(errors_are_sd){
      stat_cov <- diag((au[[i]]$Stat_err)^2)
    }else{
      stat_cov <- diag((au[[i]]$Stat_err/1.96)^2)
    }
    
    
    #Scaling pT
    pT_scaled[[i]] <- (au[[i]]$pT - min(au[[i]]$pT))/(max(au[[i]]$pT) - min(au[[i]]$pT))
    if(est_err_ell){
      err_mod <- rgasp(pT_scaled[[i]],au[[i]]$RAA_exp, kernel_type = "pow_exp",alpha = 1.9)
      ell_err <- 1/err_mod@beta_hat
    }else if(is.null(ell_err)){
      stop('Must pick a value for ell if you do\'nt estimate it')
    }
    #
    
    print(ell_err)
    if(errors_are_sd){
      sys_cov <- outer(au[[i]]$Sys_err,au[[i]]$Sys_err)*
        cov_mat(x = pT_scaled[[i]],ell = ell_err, lambda = 1,alpha = 1.9)
    }else{
      sys_cov <- outer(au[[i]]$Sys_err/1.96,au[[i]]$Sys_err/1.96)*
        cov_mat(x = pT_scaled[[i]],ell = ell_err, lambda = 1,alpha = 1.9)
    }
    
    
    covs[[i]] <- stat_cov + sys_cov
    
    #Take the experimental data, train a GP to it
    #Find the MLE length parameter when the standard deviations are fixed
    
    cur_concat <- dcast(mod_dat,variable~pT)
    colnames(cur_concat) = c('design',paste0(current_experiment,'_',colnames(cur_concat)[-1]))
    
    if(write_csvs){
      #Write csv of each dataset for easily manipulation later
      write.csv(cur_concat[,-1],file=paste0(folder,current_dset,".csv"),row.names = FALSE)
    }    
    if(first_concat){
      first_concat = FALSE
      all_mod_dat = cur_concat
      all_exp_dat = exp_dat
    }else{
      all_mod_dat = cbind(all_mod_dat,cur_concat[,-1])
      all_exp_dat = c(all_exp_dat,exp_dat)
    }
  }
  
  
  block_covs = as.matrix(bdiag(covs))
  
  
  return(list(exp_dat = all_exp_dat,
              mod_dat = all_mod_dat[,-1],
              scaled_d = scaled_d,
              block_covs = block_covs,
              ranges = ranges,
              original_dset_list = au,
              pT_scaled = pT_scaled,
              ell_err = ell_err,
              design = design))
  
}



#Covariance Function - Currently Squared Exponential
cov_mat <- function(x,ell,lambda,alpha = 2, nugget=0.){
  #inds <- 1:length(x)
  
  out_mat <- lambda^(-1)*(exp(-as.matrix(dist(x)/ell)^alpha) + 
                            nugget*diag(length(x)))
  
  return(out_mat)
}

# cov_mat <- function(pred_pts,
#                     design_pts,
#                     ell_vec,
#                     lambda,
#                     alpha_vec = rep(1.9,dim(design_pts)[2]), nugget=0.){
# 
#   M = dim(design_pts)[1]
#   N = dim(pred_pts)[1]
# 
#   if(dim(pred_pts)[2]!=dim(design_pts)[2]) stop('You got dimension problems, big fella')
# 
#   K = dim(pred_pts)[2]
# 
#   #Old one
#   # out_cov = exp(-as.matrix(1/ell_vec[1]*abs(outer(pred_pts[,1],design_pts[,1],FUN = "-")))^alpha_vec[1])
#   # if(K>1){
#   #   for(k in 2:K){
#   #     out_cov = out_cov*exp(-as.matrix(1/ell_vec[k]*abs(outer(pred_pts[,k],design_pts[,k],FUN = "-")))^alpha_vec[k])
#   #   }
#   # }
# 
#   out_cov = Exp.cov(pred_pts[,1],design_pts[,1],theta = ell_vec[1],p = alpha_vec[1])
#   if(K>1){
#     for(k in 2:K){
#       out_cov = out_cov*Exp.cov(pred_pts[,k],design_pts[,k],theta = ell_vec[k],p = alpha_vec[k])
#     }
#   }
# 
#   if(dim(pred_pts)[1] == dim(design_pts)[1]){
#     out_cov = out_cov + nugget*diag(M)
#   }
#   out_cov = (1/lambda)*out_cov
#   return(out_cov)
# 
# }





#Covariance Function - Currently Squared Exponential
#Same as commented-out cov_mat above
cov_mat_multi <- function(pred_pts,
                          design_pts,
                          ell_vec,
                          lambda,
                          alpha_vec = rep(1.9,dim(design_pts)[2]), nugget=0.){
  
  (M = dim(design_pts)[1])
  (N = dim(pred_pts)[1])
  
  if(dim(design_pts)[2]!=dim(pred_pts)[2]) stop('You got dimension problems, big fella')
  
  (K = dim(pred_pts)[2])
  
  (sig_12 = exp(-as.matrix(1/ell_vec[1]*abs(outer(pred_pts[,1],design_pts[,1],FUN = "-")))^alpha_vec[1]))
  if(K>1){
    for(k in 2:K){
      (sig_12 = sig_12*exp(-as.matrix(1/ell_vec[k]*abs(outer(pred_pts[,k],design_pts[,k],FUN = "-")))^alpha_vec[k]))
    }
  }
  
  if(dim(pred_pts)[1] == dim(design_pts)[1]){
    sig_12 = sig_12 + nugget*diag(M)
  }
  sig_12 = (1/lambda)*sig_12
  return(sig_12)
  
}

calc_mu_st_concat <- function(theta_st,
                              design,
                              model_list,
                              train_z){
  out_vec = do.call(c,lapply(seq_along(1:length(model_list)),function(i){
    comp_mod = model_list[[i]]
    z_c = train_z[[i]]
    cor_12 = cov_mat_multi(theta_st,
                           design,
                           ell_vec = 1/comp_mod@beta_hat,
                           lambda = 1,
                           alpha_vec = comp_mod@alpha,
                           nugget = comp_mod@nugget)
    #sig_22_inv = 1/comp_mod@sigma2_hat*chol2inv(chol(comp_mod@L%*%t(comp_mod@L)))
    R_inv = chol2inv(chol(comp_mod@L%*%t(comp_mod@L)))
    return(comp_mod@theta_hat + cor_12%*%R_inv%*%(z_c - comp_mod@theta_hat))
  }))
  return(out_vec)
}

calc_sig_st_concat <- function(theta_st,
                               design,
                               model_list){
  
  out_vec = do.call(c, lapply(model_list,function(comp_mod){
    #sig_11 = comp_mod@sigma2_hat
    cor_12 = cov_mat_multi(theta_st,
                           design,
                           ell_vec = 1/comp_mod@beta_hat,
                           lambda = 1,
                           alpha_vec = comp_mod@alpha,
                           nugget = comp_mod@nugget)
    L_inv =  backsolve(r = comp_mod@L, x = diag(ncol(comp_mod@L)),
                       upper.tri = FALSE)
    
    #sig_22_inv =  1/comp_mod@sigma2_hat*t(L_inv)%*%L_inv
    R_inv = t(L_inv)%*%L_inv
    
    h_st = 1
    h_d = t(t(rep(1,dim(design)[1])))
    
    n = dim(design)[1]
    q = dim(design)[2]
    (my_sig2hat <- 1/(n-1)*(Z[[1]] - comp_mod@theta_hat)%*%R_inv%*%(Z[[1]] - comp_mod@theta_hat))
    
    (marg_extra <- t(h_st - t(h_d)%*%R_inv%*%t(cor_12))%*%
        solve(t(h_d)%*%R_inv%*%h_d)%*%
        (h_st - t(h_d)%*%R_inv%*%t(cor_12)))
    #sig_22_inv = 1/comp_mod@sigma2_hat*chol2inv(chol(comp_mod@L%*%t(comp_mod@L)))
    return(comp_mod@sigma2_hat*(1 - tcrossprod(cor_12%*%R_inv,cor_12) + 
                                  marg_extra))
  })
  )
  return(out_vec)
}

