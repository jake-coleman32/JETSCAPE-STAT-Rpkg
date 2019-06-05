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

col_alpha <- function(col_str, alpha){
  if(alpha<0|alpha>1){
    stop('Alpha must be within (0,1)')
  }
  rgb_mat <- col2rgb(col_str,alpha = TRUE)
  return(rgb(rgb_mat[1],rgb_mat[2],rgb_mat[3],alpha = 255*alpha, maxColorValue = 255))
}
