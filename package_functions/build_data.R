#### Function to take file and make datasets ####
### Also makes experimental covariance, and experimental data ###
build_data <- function(folder,
                       dsets,
                       design_file,
                       range_file = NULL,
                       write_csvs = FALSE,
                       design_names = NULL,
                       est_err_ell = TRUE,
                       ell_err = NULL,
                       subset_high_pT_pbpb = FALSE,
                       errors_are_sd = TRUE,
                       add_header = FALSE,
                       split_pb_error = FALSE,
                       full_cor_au = FALSE){
  
  
  #Covariance Function - Currently Squared Exponential
  cov_mat <- function(x,ell,lambda,alpha = 2, nugget=0.){
    #inds <- 1:length(x)
    
    out_mat <- lambda^(-1)*(exp(-as.matrix(dist(x)/ell)^alpha) + 
                              nugget*diag(length(x)))
    
    return(out_mat)
  }
  
  
  #Loading design input
  #Assume the first line is the column names
  (design <- read.table(design_file,header = TRUE))
  print(paste('Colnames of design input are',paste0(colnames(design),collapse = ',')))
  
  #Dimensions of design input
  P <- dim(design)[2]
  m <- dim(design)[1]
  
  
  #Loading or creating ranges for design input
  if(!is.null(range_file)){
    load(range_file)
  }else{
    print('Note: building ranges from design')
    ranges <- vector('list', P)
    for(p in 1:P){
      ranges[[p]] <- c(min(design[,p]),max(design[,p]))
    }
  }
  
  #Design input scaled to lie in [0,1]^P
  scaled_d <- matrix(0, m, P) %>%
    as.data.frame()
  for(p in 1:dim(design)[2]){
    scaled_d[,p] <- (design[,p] - ranges[[p]][1])/(ranges[[p]][2] - ranges[[p]][1])
  }
  
  #Setting names for ranges and scaled_d
  if(!is.null(design_names)){
    print(paste0('Using design_names: ',paste0(design_names,collapse = ', ')))
    print('Note: these are also the names for ranges')
  }else{
    design_names <- paste0(colnames(design),collapse = ',')
    print('Warning: Using colnames of design_file for scaled_d and ranges')
    print(paste0('Using design_names: ', design_names))
  }
  colnames(scaled_d) <- design_names
  names(ranges) <- design_names
  
  #################################################
  ### Build covariance matrix and design output ###
  #################################################
  
  covs <- output_dset <- pT_scaled <- vector('list',length(dsets))
  
  first_concat = TRUE
  for(k in 1:length(dsets)){
    (current_dset = dsets[k])
    (current_experiment = strsplit(current_dset,'-')[[1]][1])
    
    is_PbPb = str_detect(current_experiment,'PbPb')
    is_AuAu = str_detect(current_experiment,'AuAu')
    if(!(is_PbPb + is_AuAu)){
      stop('Something wrong, neither AuAu nor PbPb')
    }else if(is_PbPb & is_AuAu){
      stop('Something wrong, both AuAu and PbPb')
    }
    
    output_dset[[k]] <- read.table(paste0(folder,current_dset,".dat"),header = add_header)
    (names(output_dset[[k]]) <- c("pT","RAA_exp","Stat_err","Sys_err",paste0("RAA_",as.character(1:m))))
    
    #Only for MATTER, where we subset to PbPb for pT > 30
    if(subset_high_pT_pbpb & is_PbPb){
      print(paste0(current_dset,', subsetting high pT'))
      output_dset[[k]] <- output_dset[[k]][which(output_dset[[k]]$pT>30),]
    }
    
    
    if(write_csvs){
      #Write csv of experimental data for easily manipulation later
      write.csv(output_dset[[k]][,1:4],file=paste0(folder,current_dset,"_exp.csv"),row.names = FALSE)
    }
    
    #Separate output, change from wide to long
    mod_dat <- dplyr::select(output_dset[[k]], -c(RAA_exp,Stat_err,Sys_err)) %>%
      melt(id = "pT") %>%
      arrange(by = pT)
    
    (exp_dat <- output_dset[[k]]$RAA_exp)
    
    
    if(errors_are_sd){
      stat_cov <- diag((output_dset[[k]]$Stat_err)^2)
    }else{
      stat_cov <- diag((output_dset[[k]]$Stat_err/1.96)^2)
    }
    
    
    #Scaling pT
    pT_scaled[[k]] <- (output_dset[[k]]$pT - min(output_dset[[k]]$pT))/
      (max(output_dset[[k]]$pT) - min(output_dset[[k]]$pT))
    
    
    #### Building Covariance Matrix ####
    if(est_err_ell){ ### This is probably defunct now
      err_mod <- rgasp(pT_scaled[[k]],output_dset[[k]]$RAA_exp, kernel_type = "pow_exp",alpha = 1.9)
      ell_err <- 1/err_mod@beta_hat
    }else if(is.null(ell_err)){
      stop('Must pick a value for ell if you do\'nt estimate it')
    }
    
    if(errors_are_sd){
      sigmas <- output_dset[[k]]$Sys_err
    }else{
      sigmas <- output_dset[[k]]$Sys_err/1.96
    }
    
    ## Building covariance matrix from given ell
    print(ell_err)
    
    ##AuAu systematic correlation
    if(is_AuAu){
      if(full_cor_au){
        sys_cov <- outer(sigmas, sigmas)
      }else{
        sys_cov <- outer(sigmas, sigmas)*
          cov_mat(x = pT_scaled[[k]],ell = ell_err, lambda = 1,alpha = 1.9)
      }
    }else if(is_PbPb){ #for PbPb, decide to split the error or not
      if(split_pb_error){
        sys_cov <- 0.5*outer(sigmas, sigmas)*
          cov_mat(x = pT_scaled[[k]],ell = ell_err, lambda = 1,alpha = 1.9) +
          
          0.5*outer(sigmas, sigmas)
        
      }else{
        sys_cov <- outer(sigmas, sigmas)*
          cov_mat(x = pT_scaled[[k]],ell = ell_err, lambda = 1,alpha = 1.9)
      }
    }else{
      stop('Problem, neither AuAu nor PbPb')
    }
    
    
    covs[[k]] <- stat_cov + sys_cov
    
    ## Current concatenated output dataset
    cur_concat <- dcast(mod_dat,variable~pT)
    colnames(cur_concat) = c('design',paste0(current_experiment,'_',colnames(cur_concat)[-1]))
    
    if(write_csvs){
      #Write csv of each dataset for easily manipulation later
      write.csv(cur_concat[,-1],file=paste0(folder,current_dset,".csv"),row.names = FALSE)
    }    
    
    ## If first dataset, create all_mod_dat,
    ### else, concatenate with existing data
    if(first_concat){
      first_concat = FALSE
      all_mod_dat = cur_concat
      all_exp_dat = exp_dat
    }else{
      all_mod_dat = cbind(all_mod_dat,cur_concat[,-1])
      all_exp_dat = c(all_exp_dat,exp_dat)
    }
  }# End k loop
  
  
  block_covs = as.matrix(bdiag(covs))
  
  
  return(list(y_exp = all_exp_dat,
              Y = all_mod_dat[,-1],
              scaled_design_input = scaled_d,
              Sigma_E = block_covs,
              ranges = ranges,
              original_dset_list = output_dset,
              pT_scaled = pT_scaled,
              ell_err = ell_err,
              design_input = design))
  
}
