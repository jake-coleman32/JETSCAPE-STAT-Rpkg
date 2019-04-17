library(MASS)
library(cluster)
library(emdbook)
library(dplyr)
library(reshape2)
library(RobustGaSP)
library(mvtnorm)
library(Matrix)
library(stringr)
library(ggplot2)
library(GGally)
library(MCMCglmm)

make_global <- function(list_obj){
  list2env(list_obj,globalenv())
}
plot_text <- function(text_to_plot){
  old_par <- par()
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, text_to_plot,
       cex = 1.6, col = "black")
  par(mar = old_par$mar)
}
col_alpha <- function(col_str, alpha){
  if(alpha<0|alpha>1){
    stop('Alpha must be within (0,1)')
  }
  rgb_mat <- col2rgb(col_str,alpha = TRUE)
  return(rgb(rgb_mat[1],rgb_mat[2],rgb_mat[3],alpha = 255*alpha, maxColorValue = 255))
}

setwd("/Users/Jake/Box/Research/Computer_Emulation/Shanshan_Project/")

source("/Users/Jake/Box/Research/Computer_Emulation/Shanshan_Project/Code/prepare_concat_data.R")
source('/Users/Jake/Box/Research/Computer_Emulation/Shanshan_Project/Code/assess_concat_analysis.R')
source('/Users/Jake/Box/Research/Computer_Emulation/Shanshan_Project/Code/training_prediction_concat.R')


folder = "Data/to_Jake/MATTER/forSTAT-MATTER/"
system(paste0('ls ',folder))
design_file = paste0(folder,'matter_lhc.txt')
range_file = paste0(folder,'ranges.Rdata')

#either MATTER/LBT/MATLBT1/MATLBT1_separate_tune
#RAA_prefix = 'MATLBT1_ell_10'
#RAA_folder = 'Results/RAA_from_median_paras_ell_p2/'

all_dsets <- c(
  "AuAu200-cen-00-10"
  ,"AuAu200-cen-40-50"
  ,"PbPb2760-cen-00-05"
  ,"PbPb2760-cen-30-40"
  ,"PbPb5020-cen-00-10"
  ,"PbPb5020-cen-30-50"
)

RAA_dsets <- c(
  "AuAu200-cen-00-10"
  ,"AuAu200-cen-40-50"
  ,"PbPb2760-cen-00-05"
  ,"PbPb2760-cen-30-40"
  ,"PbPb5020-cen-00-10"
  ,"PbPb5020-cen-30-50"
)

all_systems = c("AuAu200"
                ,"PbPb2760"
                ,"PbPb5020"
)
all_data_names = c('PHENIX',
                   'ATLAS',
                   'CMS')

elt = 'all'

if(elt=='all'){
  sys_index = 1:3
  collider = 'Combined'
}else if(elt == 'au'){
  sys_index =  1
  collider = 'RHIC'
}else if (elt == 'pb'){
  sys_index = 2:3
  collider = 'LHC'
}

systems_to_calibrate = all_systems[sys_index]
data_names = all_data_names[sys_index]



if(length(systems_to_calibrate)==3){
  (hist_main = "All Datasets Simultaneous Calibration")
}else{
  (hist_main = paste(paste0(systems_to_calibrate,collapse = " & "),"Calibration"))
} 
(dsets_to_use = which(str_detect(all_dsets, paste(systems_to_calibrate,collapse = "|"))))
(dset_strings = all_dsets[dsets_to_use])
(RAA_med_strings <- RAA_dsets[dsets_to_use])

(renumbered_dsets_to_use = which(str_detect(dset_strings, paste(systems_to_calibrate,collapse = "|"))))

#### Some String things###
collision_sys_names = systems_to_calibrate
(cen_groups <- split(renumbered_dsets_to_use,rep(systems_to_calibrate,each=2)))

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

(cen_names <- make_cen_names(cen_groups,  dset_strings))

(cen_groups_design <- setNames(as.list(renumbered_dsets_to_use), rep(systems_to_calibrate,each=2)))
(cen_names_design <- make_cen_names(cen_groups_design, dset_strings))

make_RAA_med_vals <- function(cen_groups, RAA_med_strings,
                              RAA_folder, RAA_prefix, cen_names){
  med_vals <- vector('list',length(cen_groups))
  (names(med_vals) <- names(cen_groups))
  
  for(j in 1:length(cen_groups)){
    med_vals[[j]] <- vector('list',length(cen_groups[[j]]))
    (names(med_vals[[j]]) <- cen_names[[j]])
    
    for(k in 1:length(cen_groups[[j]])){
      
      #R_AA Median Values
      print(paste0(RAA_prefix,'_',RAA_med_strings[cen_groups[[j]][k]]))
      (pt_RAA <- read.table(paste0(RAA_folder,RAA_prefix,'_',RAA_med_strings[cen_groups[[j]][k]],'.dat')))
      
      #Subset to high pT for MATTER
      indices_to_use = 1:dim(pt_RAA)[1]
      if(RAA_prefix=='MATTER'){
        if(names(med_vals)[j]!="AuAu200"){
          indices_to_use = which(pt_RAA[,1]>30)
        }
      }
      med_vals[[j]][[k]] <- pt_RAA[indices_to_use,2]
      
    }
  }
  
  print('MAKE SURE: folder is')
  print(RAA_folder)
  print('prefix is')
  print(RAA_prefix)
  
  return(med_vals)
  
}


med_vals <- make_RAA_med_vals(cen_groups,RAA_med_strings, RAA_folder,
                              RAA_prefix, cen_names)


med_vals_length <- do.call(c,lapply(med_vals,function(x){do.call(c,lapply(x,length))}))
dset_length <- do.call(c,lapply(comp_mod$all_data$original_dset_list,function(x){length(x$pT)}))

if(sum(med_vals_length != dset_length)){
  print('NEED TO DO STUPID SUBSETTING')
}

##STUPID SUBSETTING
if(RAA_prefix=='MATTER'){
  for(i in 1:2){
    med_vals[[1]][[i]] <- med_vals[[1]][[i]][-1]
  }
}else{
  med_vals_chop <- lapply(med_vals,function(sys){
    lapply(sys,function(cen){
      cen[-1]
    })
  })
  med_vals <- med_vals_chop
}



q_run <- 3
kernel_type_run <- 'matern_5_2'
nugget_est_run <- TRUE
nugget_run = 0
subset_for_MATTER = TRUE
ell_run = 10
est_err_ell_run = FALSE
header_for_LBT = FALSE

if(subset_for_MATTER){
  print('WARNING: SUBSETTING FOR MATTER')
}
if(header_for_LBT){
  print('WARNING: INCLUDING HEADER, DONT LOSE FIRST DATA POINT')
}


comp_mod = train_models(dset_path = folder,
                        datasets = all_dsets[dsets_to_use],
                        design_file = design_file,
                        range_file = range_file,
                        write_csvs = FALSE,
                        est_err_ell = est_err_ell_run,
                        ell_err = ell_run,
                        subset_high_pT_pbpb = subset_for_MATTER,
                        errors_are_sd = TRUE,
                        add_header = header_for_LBT,

                        design_names = c('A+C','A/(A+C)', 'D', 'Q'),
                        #design_names = c('A','B','C','D'),
                        q = q_run,
                        
                        nugget = nugget_run,
                        nugget.est = nugget_est_run,
                        kernel_type = kernel_type_run
)

#nuggets
(nug_est <- as.character(formatC(do.call(c,lapply(comp_mod$train_mods,function(x)x@nugget)),format = "e", digits = 2)))


Sigma_E <- comp_mod$all_data$block_covs
y_E_concat <- comp_mod$all_data$exp_dat
(D_c_theta <- comp_mod$all_data$scaled_d)
(design <- comp_mod$all_data$design)
ranges<- comp_mod$all_data$ranges
ell_err <- comp_mod$all_data$ell_err
original_dset_list = comp_mod$all_data$original_dset_list
J <- length(original_dset_list)


save_list <- list()
save_list$comp_mod <- comp_mod
save_list$dset_strings <- dset_strings



#Rotate with PCA
Y <- comp_mod$all_data$mod_dat
y_sd = comp_mod$y_sd
y_means = comp_mod$y_means


Y_final <- comp_mod$Y_final

#Some exploratory measures from assess_concat_analysis.R
normality_test(Y_final)
(Fq <- calc_Vq(comp_mod$Y_svd,plotit = TRUE))[1:6]

q <- length(comp_mod$train_mods)
V <-  comp_mod$Y_svd$V[,1:q]

save_list$V = V

#############
##Prediction
############

save_list$final_mod = comp_mod$final_mod
# comp_mod = final_mod[[1]]

#
#############
#Calibration
#############


mh_cal <- function(niters = 1E4,
                   burnin_prop = 0.3,
                   mod = comp_mod,
                   t_kap = 0.1,
                   proposal_cor_mat =diag(dim(D_c_theta)[2]),
                   upper_theta_trunc = 1,
                   bad_corners = list()){#Keep this for now, but never turn it on
  
  #Do various sampler setups
  time_start <- proc.time()
  
  burnin <- burnin_prop*niters
  
  t_out <- matrix(0,niters,dim(proposal_cor_mat)[2])
  t_ratio <- ob_ratio <- 0
  
  #draw from priors
  (t_cur <- t(runif(dim(proposal_cor_mat)[2],0,upper_theta_trunc)))
  
  #Get current values
  (X_cur <- t_cur)
  colnames(X_cur) <- colnames(D_c_theta)
  cur_preds = exp_predictions(theta = X_cur,
                              model_list = mod$train_mods,
                              y_sd = mod$y_sd,
                              y_means = mod$y_means,
                              Y_svd = mod$Y_svd,
                              q = length(mod$train_mods),
                              cov_extra = mod$cov_extra)
  
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
    
    
    X_st <- t_st
    colnames(X_st) <- colnames(D_c_theta)
    
    
    st_preds = exp_predictions(theta = X_st,
                               model_list = mod$train_mods,
                               y_sd = mod$y_sd,
                               y_means = mod$y_means,
                               Y_svd = mod$Y_svd,
                               q = length(mod$train_mods),
                               cov_extra = mod$cov_extra)
    
    y_st_means <- st_preds$mean
    y_st_cov <- st_preds$cov
    
    y_cur_means <- cur_preds$mean
    y_cur_cov <- cur_preds$cov
    
    
    #Get likelihoods
    
    
    (like_st <- mvtnorm::dmvnorm(x = as.numeric(y_E_concat),
                                 mean = as.numeric(y_st_means),
                                 sigma = Sigma_E + y_st_cov, #cov_extra_phys_cal is part of y_st_cov, and y_cur_cov
                                 log = TRUE))
    (like_cur <- mvtnorm::dmvnorm(x = as.numeric(y_E_concat),
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
  
  return(list(params = t_out))
  
}

res <- mh_cal(mod = comp_mod, niters = 1E5,t_kap = 7E-3
              #, upper_theta_trunc = 1
               # ,bad_corners = list(list(c('high',NA,'low',NA),
               #                         c(0.82,NA,0.2,NA)),
               #                     list(c('low',NA,'high',NA),
               #                          c(0.33,NA,0.9,NA)))
)

save_list$res <- res
#proposal_cor_mat = matrix(c(1,-.8,-.8,1),ncol = 2))
#save(res,file = paste0(save_path,'res_old.Rdata'))

plot_text(paste0(kernel_type_run,', ',q,' PCs, ',nugget_est_run, ' nugget est, \n',
                 formatC(dim(res$params)[1], format = "e", digits = 1),' iters, \n',
                 hist_main))

param_plot <- matrix(0,dim(res$params)[1],dim(res$params)[2])
#param_plot <- matrix(0,dim(res$params)[1],dim(res$params)[2])
for(j in 1:dim(param_plot)[2]){
  param_plot[,j] <- res$params[,j]*(ranges[[j]][2] - ranges[[j]][1]) + ranges[[j]][1]
}
#param_plot_sqrt <- sqrt(param_plot)
#names(ranges) <- c("A","B","C","D")
(colnames(param_plot) <- names(ranges))


(save_suffix = paste0("_matter_pb_ell_10_", kernel_type_run,'_',q))
save_list$save_suffix <- save_suffix

save_list$ranges = ranges
save_list$param_plot = param_plot


save_path <- "Results/long_ell_4_4/MATTER/"

plot_heatmaps(param_plot,ranges,pairs = list(c(1,2),c(3,4)),
              param_names = colnames(param_plot),
              save_pics = TRUE,
              save_param_names = list('A_C','B_D'),
              save_path = save_path,
              save_end = save_suffix,
              plt_points = TRUE,design = design)


# save_pics = FALSE
# if(save_pics) pdf(paste0(save_path,'Q_hist',save_suffix,'.pdf'))
hist(param_plot[,5], xlab = 'Q',main = 'Posterior histogram of Q')
# if(save_pics) dev.off()



###Single heat_pairplot
# heat_pairplot(param_plot,
#               ranges,
#               split_num = dim(param_plot)[1],
#               param_names = colnames(param_plot), 
#               title = 'Calibration Pairplot',
#               save_pics = FALSE,
#               save_path = save_path,
#               save_end = save_suffix,
#               include_legend = FALSE,
#               legend_names = NULL,
#               legend_cols = NULL,
#               off_diag = "heatmap")


transforming = TRUE
plotting = FALSE
if(transforming){
  transformed_params <- transform_draws(param_plot = param_plot,
                                        cols_to_transform = c(1,2),
                                        new_names = c('A','C'))
  
  ranges_t <- transform_ranges(ranges = ranges,
                               vars_to_transform = c(1,2),
                               new_names = c('A','C'))
  
  if(plotting){
    plot_heatmaps(transformed_params,ranges_t,pairs = list(c(1,2)),
                  param_names = colnames(transformed_params),
                  save_pics = TRUE,
                  save_param_names = list('A_C'),
                  save_path = save_path,save_end = paste0(save_suffix,'_transformed'),
                  plt_points = FALSE,design = design)
  }
  
  save_list$transformed_params <- transformed_params
  save_list$ranges_t <- ranges_t
}

save_path
save_suffix

save(save_list,file=paste0(save_path,'save_list',save_suffix,'.Rdata'))

model_specs <- paste0(hist_main,' with q = ',q,', ',kernel_type_run,' kernel, ', nugget_est_run, ' nugget est, ',
                      formatC(dim(res$params)[1], format = "e", digits = 1),' iters')
(model_specs <- c(hist_main, 
                  as.character(q),
                  kernel_type_run,
                  ell_run,
                  nugget_est_run,
                  as.character(formatC(dim(res$params)[1], format = "e", digits = 1)),
                  paste0(nug_est,collapse = ' '),
                  folder,
                  as.character(Sys.Date())))
spec_types <- c('Dsets','PCs','Ell', 'Kernel','nugget_est','Iterations','Nuggets','Data Folder','Date')

specs_out <- paste(spec_types,model_specs,sep = ': ',collapse = '\n\n')

#BE CAREFUL SAVE PATH IS CORRECT
save_path
specs_end = '_matter_pb_ell_10'
write(as.character(specs_out),file = paste0(save_path,paste0('model_specs',specs_end,'.txt')))



#####
#Extra visualization/analysis
#####

# setwd("/Users/Jake/Box/Research/Computer_Emulation/Shanshan_Project/")
# save_path <-"Results/LBT+MATTER_correction/LBT+MATTER-Method-1/pow_exp_5_PCs/"
# setwd(save_path)
# load('save_list_m1_pow_exp_5.Rdata')
# make_global(save_list)
# setwd("/Users/Jake/Box/Research/Computer_Emulation/Shanshan_Project/")


J <- length(comp_mod$all_data$original_dset_list)
delta_means = lapply(seq_along(1:J),function(x){0})

## Just for testing ##
# load(paste0(save_path,'save_list_lbt_pow_exp_all.Rdata'))
# make_global(save_list)
# res <- list()
# res$params <- unmake_param_plot(param_plot,ranges)
########

#Posteriors
plot_draws_together(draws = res$params,
                    train_mod = comp_mod$train_mods,
                    rot_mat = comp_mod$Y_svd$v[,1:comp_mod$q],
                    original_dset_list = comp_mod$all_data$original_dset_list,
                    scale_sd_vec = comp_mod$y_sd,
                    scale_mean_vec = comp_mod$y_means,
                    dset_labels = collision_sys_names,
                    data_orig_labels = data_names,
                    cen_groups = cen_groups,
                    cen_labels = cen_names,
                    
                    q_index = 5,
                    condition_on_q = FALSE,
                    
                    med_pred_type = 'emulator', ##RAA or emulator
                    RAA_med_vals = med_vals,
                    
                    title_end = 'Posterior',
                    # col_str_vec = c('skyblue'),
                    # col_pts_vec = c('darkblue'),
                    
                    col_str_vec = c('skyblue','pink'),
                    col_pts_vec = c('darkblue','darkred'),
                    
                    
                    alpha_val = 0.05,
                    save_pics = TRUE,
                    save_path = save_path,
                    delta_list = delta_means,
                    #save_end = paste0(save_suffix,"_all_pts_posterior"),
                    save_end = paste0(save_suffix,"_posterior"),
                    legend_spots = replicate(J,'topleft'),
                    include_med = FALSE,
                    include_arrows = TRUE)

#Design
plot_draws_together(draws = comp_mod$all_data$scaled_d[1:60,],
                    num_samples = dim(comp_mod$all_data$scaled_d)[1],
                    train_mod = comp_mod$train_mods,
                    rot_mat = comp_mod$Y_svd$v[,1:comp_mod$q],
                    original_dset_list = comp_mod$all_data$original_dset_list,
                    scale_sd_vec = comp_mod$y_sd,
                    scale_mean_vec = comp_mod$y_means,
                    
                    ###Which one? Together or separate?
                    
                    #Separate
                    #alpha_val = 1,
                    # dset_labels = paste0(names(cen_groups_design),paste0("_",c(1:2))),##Hacky - would override dset_labels
                    # cen_groups = cen_groups_design,
                    # cen_labels = cen_names_design,
                    # col_str_vec = c('pink'),
                    # col_pts_vec = c('darkred'),
                    
                    #Together
                    alpha_val = 0.8,
                    dset_labels = collision_sys_names,
                    cen_groups = cen_groups,
                    cen_labels = cen_names,
                    col_str_vec = c('skyblue','pink'),
                    col_pts_vec = c('darkblue','darkred'),
                    
                    title_end = 'Design',
                    
                    
                    save_pics = FALSE,
                    save_path = save_path,
                    delta_list = delta_means,
                    #save_end = paste0(save_suffix,"_all_pts_posterior"),
                    save_end = paste0(save_suffix,"_design"),
                    legend_spots = replicate(J,'topleft'),
                    include_legend = TRUE,
                    include_med = FALSE,
                    include_arrows = TRUE)

# plot_draws(draws = res$params,
#            train_mod = comp_mod$train_mods,
#            rot_mat = comp_mod$Y_svd$v[,1:comp_mod$q],
#            original_dset_list = comp_mod$all_data$original_dset_list,
#            scale_sd_vec = comp_mod$y_sd,
#            scale_mean_vec = comp_mod$y_means,
#            title_end = 'Posterior',
#            alpha_val = 0.03,
#            save_pics = FALSE,
#            save_path = save_path,
#            delta_list = delta_means,
#            #save_end = paste0(save_suffix,"_all_pts_posterior"),
#            save_end = paste0(save_suffix,"_posterior"),
#            legend_spots = replicate(J,'topleft'),
#            include_med = TRUE)
# plot_draws(draws = matrix(runif(dim(res$params)[1]*dim(res$params)[2]),
#                           ncol = dim(res$params)[2]),
#            train_mod = comp_mod$train_mods,
#            rot_mat = comp_mod$Y_svd$v[,1:comp_mod$q],
#            original_dset_list = comp_mod$all_data$original_dset_list,
#            scale_sd_vec = comp_mod$y_sd,
#            scale_mean_vec = comp_mod$y_means,
#            dset_labels = collision_sys_names,
#            title_end = 'Prior',
#            alpha_val = 0.03,
#            save_pics = TRUE,
#            save_path = save_path,
#            delta_list = delta_means,
#            #save_end = paste0(save_suffix,"_all_pts_prior"),
#            save_end = paste0(save_suffix,"_prior"),
#            legend_spots = replicate(J,'topleft'))
# 
# plot_draws(draws = comp_mod$all_data$scaled_d,
#            train_mod = comp_mod$train_mods,
#            rot_mat = comp_mod$Y_svd$v[,1:comp_mod$q],
#            original_dset_list = comp_mod$all_data$original_dset_list,
#            scale_sd_vec = comp_mod$y_sd,
#            scale_mean_vec = comp_mod$y_means,
#            dset_labels = collision_sys_names,
#            title_end = "Design",
#            alpha_val = .8,
#            num_samples = dim(comp_mod$all_data$scaled_d)[1],
#            save_pics = FALSE,
#            save_path = save_path,
#            save_end = "design",
#            delta_list = delta_means,
#            include_legend = FALSE)

##### Extra for Shanshan
q_index = 5
q_mode <- density(transformed_params[,q_index])$x[which.max(density(transformed_params[,q_index])$y)]
param_to_summarize = transformed_params

condition_on_q = FALSE
if(condition_on_q){
  q_width = 0.3
  good_q <- which(param_to_summarize[,q_index] > (q_mode - q_width) &
                    param_to_summarize[,q_index] < (q_mode + q_width) )
  param_to_summarize <- param_to_summarize[good_q, ]
}

rand_sample = sample(1:dim(param_to_summarize)[1],1E4)
write.table(round(param_to_summarize[rand_sample,],3),
            file = paste0(save_path,'calibration_draws',save_suffix,'.txt'),
            row.names = FALSE)

(summary_stats = apply(param_to_summarize,2,function(x){
  med = median(x)
  dens_x <- density(x)$x
  dens_y <- density(x)$y
  (mode = dens_x[which.max(dens_y)])
  mean = mean(x)
  return(round(c(median = med, mean = mean,mode = mode),3))
}))

(q_summary <- round(quantile(param_to_summarize[,q_index],probs = c(0.05,0.5,0.95)),2))

write.table(summary_stats,
            file = paste0(save_path,'summary_stats',save_suffix,'.txt'))


#######################
### Making combined pairplot
#########################
save_path <- "Results/long_ell_4_4/MATTER+LBT_1/"
system(paste0('ls ',save_path,'*.Rdata'))
all_suffixes <- c('m1_ell_10_matern_5_2_5','m1_au_ell_10_matern_5_2_5',
                  'm1_pb_ell_10_matern_5_2_5')
labels <- c('Combined','RHIC','LHC')

# all_suffixes <- c('lbt_ell_10_matern_5_2_3', 'matter_ell_10_matern_5_2_3')
# labels <- c('LBT','MATTER')


(all_lists <- paste0(paste0(save_path,'save_list_'),all_suffixes,'.Rdata'))

#set.seed(47)
load(all_lists[[1]])
make_global(save_list)
names(save_list)
save_suffix

make_combined_pairplot(all_lists,
                       labels,
                       ranges = ranges,
                       save_pics = TRUE,
                       save_path = save_path,
                       save_end = 'lbt_v_matter_ell_10_matern_5_2_3',
                       cols_to_use = 1:4,
                       col_names = c('A','B','C','D'),
                       scatter_alpha = 0.01,
                       num_samples = 1E4,
                       str_param_to_use = 'param_plot',
                       make_pdf = TRUE)

for(i in 1:length(dset_strings)){
  print(length(comp_mod$all_data$original_dset_list[[i]]$pT))
}
do.call(c,do.call(c,lapply(med_vals,function(x){lapply(x,function(i)length(i))})))
