### Helper functions ##

rescale_posteriors <- function(cal_posteriors,
                               ranges){
  params <- matrix(0,dim(cal_posteriors)[1],dim(cal_posteriors)[2])

    for(j in 1:dim(params)[2]){
      params[,j] <- cal_posteriors[,j]*(ranges[[j]][2] - ranges[[j]][1]) + ranges[[j]][1]
    }
  
  if(is.null(names(ranges))){
    print('Note: No names for ranges, parameters')
  }else{
    colnames(params) <- names(ranges)
  }
  
  return(params)
}

#Make list with centrality name 
make_cen_names <- function(cen_groups, dset_strings){
  (cen_names = vector('list',length(cen_groups)))
  (names(cen_names) <- names(cen_groups))
  
  for(j in 1:length(cen_groups)){
    
    (cen_names[[j]] <- character(length(cen_groups[[j]])))
    
    for(k in 1:length(cen_groups[[j]])){
      
      #Centrality Names
      (cen_index = cen_groups[[j]][k])
      (dset_st <- dset_strings[cen_index])
      (cen_index = regexpr('cen-',dset_st)[1] + 4)
      cen_strings = substring(dset_st,first = cen_index) %>%
        str_split('-') %>%
        .[[1]] %>%
        as.numeric() %>%
        paste0('%',collapse = ' - ')
      (cen_names[[j]][k] <- cen_strings)
      
    }
  }
  return(cen_names)
}


build_labels <- function(data_file_strings,
                         system_names,
                         data_names,
                         elt = 'all'){
  if(elt=='all'){
    sys_index = 1:3
    collider = 'Combined'
  }else if(elt == 'au'){
    sys_index =  1
    collider = 'RHIC'
  }else if (elt == 'pb'){
    sys_index = 2:3
    collider = 'LHC'
  }else{
    stop('elt must be all, pb or au')
  }
  
  systems_to_calibrate = system_names[sys_index]
  (dsets_to_use = which(str_detect(data_file_strings, paste(systems_to_calibrate,collapse = "|"))))
  
  dset_strings_to_use = data_file_strings[dsets_to_use] #i.e. the files shanshan gives me
  data_names_to_use = data_names[sys_index] #i.e. CMS, ATLAS, ALICE
  
  ## Need to renumber because we took a subset of 
  (renumbered_dsets_to_use <- which(str_detect(dset_strings_to_use, paste(systems_to_calibrate,collapse = "|"))))
  
  (centrality_groups <- split(renumbered_dsets_to_use,rep(systems_to_calibrate,each=2)))
  
  
  centrality_names <- make_cen_names(centrality_groups,
                                     dset_strings_to_use)
  
  return(list(centrality_groups = centrality_groups,
              centrality_names = centrality_names,
              systems_to_calibrate = systems_to_calibrate,
              data_names_to_use = data_names_to_use,
              dset_strings_to_use = dset_strings_to_use))
  
}


make_output_list <- function(folder,
                             dsets,
                             subset_high_pT_pbpb = FALSE,
                             errors_are_sd = TRUE,
                             add_header = FALSE #LBT only
){
  
  ###################################
  ### Build design output list  ###
  ###################################
  
  output_list <- vector('list',length(dsets))
  
  first_concat = TRUE
  for(k in 1:length(dsets)){
    (current_dset = dsets[k])
    (current_experiment = strsplit(current_dset,'-')[[1]][1])
    
    is_PbPb = str_detect(current_experiment,'PbPb')
    
    output_list[[k]] <- read.table(paste0(folder,current_dset,".dat"),header = add_header)
    
    m <- dim(output_list[[k]])[2] - 4
    
    (names(output_list[[k]]) <- c("pT","RAA_exp","Stat_err","Sys_err",paste0("RAA_",as.character(1:m))))
    
    #Only for MATTER, where we subset to PbPb for pT > 30
    if(subset_high_pT_pbpb & is_PbPb){
      print(paste0(current_dset,', subsetting high pT'))
      output_list[[k]] <- output_list[[k]][which(output_list[[k]]$pT>30),]
    }
    
    
    
  }
  return(output_list)
}

#Make named colors more translucent
col_alpha <- function(col_str, alpha){
  if(alpha<0|alpha>1){
    stop('Alpha must be within (0,1)')
  }
  rgb_mat <- col2rgb(col_str,alpha = TRUE)
  return(rgb(rgb_mat[1],rgb_mat[2],rgb_mat[3],alpha = 255*alpha, maxColorValue = 255))
}
