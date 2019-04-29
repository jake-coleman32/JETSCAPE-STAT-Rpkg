## Training GP models


train_GPs <- function(design_input,
                      computer_output,
                      experimental_vector,
  
                      range_file = NULL,
                      R = 5,
                      
                      nugget = 0,
                      nugget.est = F,
                      kernel_type = 'matern_5_2'){
  
  
  #Renaming for easier management
  design <- design_input
  Y <- computer_output
  
  #Dimensions of Y is (n x q)
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  
  
  #Collect mean and standard deviation of each column of our output matrix
  y_sd = apply(Y,2,sd)
  y_means = apply(Y,2,mean)
  
  #Center and scale
  Y_final <- as.matrix(sweep(Y,2,y_means)) %*% diag(1/y_sd)
  
  #Get the SVD of centered/scaled Y
  #Use only the first R components
  Y_svd <- svd(Y_final)
  V = Y_svd$v[, 1:R]
  S = diag(Y_svd$d[1:R])
  
  #Rotate with PCA
  Z <-as.matrix(Y_final)%*%V %>%
    as.data.frame()
  
  if(is.null(range_file)){
    print('Note: Using design_input to build ranges, might underestimate ranges')
    ranges <- vector('list',P)
    for(p in 1:P){
      ranges[[p]] <- c(min(design[,p]), max(design[,p]))
    }
  }else{
    load(range_file)
    #Maybe some check that the ranges file loaded properly?
  }
  
  #Scale design to use in the RGaSP package
  scaled_design <- matrix(0,dim(design)[1],dim(design)[2])
  for(p in 1:P){
    scaled_design[,p] <- (design[,p] - ranges[[p]][1])/(ranges[[p]][2] - ranges[[p]][1])
  }
  
  
  #############
  ##Prediction
  ############
  
  final_mod <- lapply(Z,rgasp,design = scaled_design,
                      nugget.est=nugget.est,
                      nugget = nugget,
                      kernel_type = kernel_type
  )
  
  
  #Account for extra variation lost. Remember Y is a (n x q) matrix
  if(R < q){
    V_b = Y_svd$v[,(R+1):q]
    S_b = diag(Y_svd$d[(R+1):q])
    cov_extra <- 1/(n-1)*diag(y_sd)%*%V_b%*%S_b^2%*%t(V_b)%*%diag(y_sd)
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
              R = R))
}