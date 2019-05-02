#Training the GP Emulators via PCA

library(MASS)
library(cluster)
library(emdbook)
library(dplyr)
library(reshape2)
library(RobustGaSP)
library(mvtnorm)
library(Matrix)
library(stringr)


train_models <- function(dset_path,
                         datasets,
                         design_file,
                         range_file,
                         write_csvs = FALSE,
                         est_err_ell = FALSE,
                         ell_err = 0.2,
                         subset_high_pT_pbpb = FALSE,
                         add_header = FALSE,
                         errors_are_sd = TRUE,
                         design_names = NULL,
                         q = 5,
                         

                         nugget = 0,
                         nugget.est = F,
                         kernel_type = 'matern_5_2'){
  
  all_data <- concatenate_dsets(folder = dset_path,
                                dsets = datasets,
                                design_file = design_file,
                                range_file = range_file,
                                write_csvs = write_csvs,
                                est_err_ell = est_err_ell,
                                design_names = design_names,
                                ell_err = ell_err,
                                subset_high_pT_pbpb = subset_high_pT_pbpb,
                                errors_are_sd = errors_are_sd,
                                add_header = add_header)
  
  
  Sigma_E <- all_data$block_covs
  y_E_concat <- all_data$exp_dat
  D_c_theta <- all_data$scaled_d
  (ranges <- all_data$ranges)
  ell_err <- all_data$ell_err
  original_dset_list = all_data$original_dset_list
  J <- length(original_dset_list)
  
  
  
  design <- matrix(0,dim(D_c_theta)[1],dim(D_c_theta)[2])
  for(i in 1:dim(D_c_theta)[2]){
    design[,i] <- D_c_theta[,i]*(ranges[[i]][2] - ranges[[i]][1]) + ranges[[i]][1]
  }
  
  
  #Rotate with PCA
  Y <- all_data$mod_dat
  y_sd = apply(Y,2,sd)
  y_means = apply(Y,2,mean)
  
  
  Y_final <- as.matrix(sweep(Y,2,y_means)) %*% diag(1/y_sd)
  Y_svd <- svd(Y_final)
  
  #Some exploratory measures from assess_concat_analysis.R
  print(normality_test(Y_final))
  print(calc_Vq(Y_svd)[1:q],plotit = FALSE)
  
  V = Y_svd$v[,1:q]
  S = diag(Y_svd$d[1:q])
  
  
  Z <-as.matrix(Y_final)%*%V %>%
    as.data.frame()
  
  
  #############
  ##Prediction
  ############
  
  final_mod <- lapply(Z,rgasp,design = D_c_theta,
                      nugget.est=nugget.est,
                      nugget = nugget,
                      kernel_type = kernel_type
  )
                      # trend = as.matrix(cbind(rep(1,m),D_c_theta)),
                      #  nugget = 1E-2,
                      #kernel_type = 'pow_exp')
  
  
  m <- dim(Y_svd$u)[1] #I think
  if(q<dim(Y_svd$v)[2]){
    V_b = Y_svd$v[,(q+1):dim(Y_svd$v)[2]]
    S_b = diag(Y_svd$d[(q+1):dim(Y_svd$v)[2]])
    cov_extra <- 1/(m-1)*diag(y_sd)%*%V_b%*%S_b^2%*%t(V_b)%*%diag(y_sd)
  }else{
    cov_extra = 0
  }
  
  
  
  return(list(train_mods = final_mod,
              y_sd = y_sd,
              y_means = y_means,
              Y_final = Y_final,
              Y_svd = Y_svd,
              all_data = all_data,
              cov_extra = cov_extra,
              q = q))
}



exp_predictions <- function(theta,
                            model_list,
                            y_sd,
                            y_means,
                            Y_svd,
                            cov_extra,
                            q){
  
  
  mod_preds <- lapply(model_list, RobustGaSP::predict, testing_input = theta)
  
  z_pred_means <- do.call(cbind,lapply(mod_preds,function(x){x$mean}))
  z_pred_vars <- do.call(cbind,lapply(mod_preds,function(x){x$sd^2}))
  
  V <- Y_svd$v[,1:q]
  y_pred_means <- (z_pred_means%*%t(V)%*%diag(y_sd)) %>%
    sweep(2,y_means,FUN = "+")
  
  if(!is.null(dim(theta)) && (dim(theta)[1]>1)){
    print("WARNING: predicting for multiple observations: returning 
          just predictive means for 5PCs, not covariance matrix in physical space")
    y_cov = z_pred_vars
  }else{
    z_pred_vars <- as.numeric(z_pred_vars)
    y_cov <- diag(y_sd)%*%V%*%diag(z_pred_vars)%*%t(V)%*%diag(y_sd) + cov_extra
  }
  
 
  return(list(mean = t(y_pred_means), cov = y_cov))
}

