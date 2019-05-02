
# Make experimental predictions at the given calibration parameter #
# This function is done #
exp_predictions <- function(theta,
                            GP_list,
                            Y_svd,
                            y_sd,
                            y_means,
                            cov_extra = 0){
  
  
  R <- length(GP_list)
  
  # Get predictions at the input
  mod_preds <- lapply(GP_list, RobustGaSP::predict, testing_input = theta)
  
  #Get the prediction means and variances at each new point
  z_pred_means <- do.call(cbind, lapply(mod_preds, function(x){x$mean}))
  z_pred_vars <- do.call(cbind, lapply(mod_preds, function(x){x$sd^2}))
  
  #Rotate the means to the physical space
  V <- Y_svd$v[,1:R]
  y_pred_means <- (z_pred_means %*% t(V) %*% diag(y_sd)) %>%
    sweep(2, y_means, FUN = "+")
  

  
  
  # Collect the predictive variances
  # If theta is a matrix with multiple rows, only get the univariate variance
  # If theta is a single multi-dimensional point, get the rotate variance matrix for calibration
  if(!is.null(dim(theta)) && (dim(theta)[1] > 1)){
    print("WARNING: predicting for multiple observations: returning 
          just predictive means for 5PCs, not covariance matrix in physical space")
    y_cov = z_pred_vars
  }else{
    z_pred_vars <- as.numeric(z_pred_vars)
    y_cov <- diag(y_sd) %*% V %*% diag(z_pred_vars) %*% t(V) %*% diag(y_sd) + cov_extra
    
  }
  
  
  return(list(mean = t(y_pred_means), cov = y_cov))
}
