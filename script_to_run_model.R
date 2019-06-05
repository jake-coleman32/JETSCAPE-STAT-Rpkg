#Script to run analysis with new functions

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


setwd("/Users/Jake/Box/Research/Computer_Emulation/Shanshan_Project/")

source("/Users/Jake/Box/Research/Computer_Emulation/Shanshan_Project/JETSCAPE_STAT_pkg/package_functions/build_data.R")
source('/Users/Jake/Box/Research/Computer_Emulation/Shanshan_Project/JETSCAPE_STAT_pkg/package_functions/calibration_sampler.R')
source('/Users/Jake/Box/Research/Computer_Emulation/Shanshan_Project/JETSCAPE_STAT_pkg/package_functions/exp_predictions.R')
source('/Users/Jake/Box/Research/Computer_Emulation/Shanshan_Project/JETSCAPE_STAT_pkg/package_functions/train_GPs.R')
source('/Users/Jake/Box/Research/Computer_Emulation/Shanshan_Project/JETSCAPE_STAT_pkg/package_functions/helper_functions.R')
source('/Users/Jake/Box/Research/Computer_Emulation/Shanshan_Project/JETSCAPE_STAT_pkg/package_functions/plotting_functions.R')
source('/Users/Jake/Box/Research/Computer_Emulation/Shanshan_Project/JETSCAPE_STAT_pkg/package_functions/calibration_pairplot.R')
source('/Users/Jake/Box/Research/Computer_Emulation/Shanshan_Project/JETSCAPE_STAT_pkg/package_functions/plot_predictions.R')



data_folder = "Data/LBT+MATTER_correction/LBT+MATTER-Method-1/forSTAT-5D-70pts/"

all_dsets <- c(
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


labels <- build_labels(data_file_strings = all_dsets,
                       system_names = all_systems,
                       data_names = all_data_names)



system(paste0('ls ',data_folder))
design_file = paste0(data_folder,'latin_hc_five.txt')
range_file = paste0(data_folder,'ranges.Rdata')

subset_for_MATTER = FALSE
ell_run = 0.2
est_err_ell_run = FALSE
header_for_LBT = FALSE

if(subset_for_MATTER){
  print('WARNING: SUBSETTING FOR MATTER')
}
if(header_for_LBT){
  print('WARNING: INCLUDING HEADER, DONT LOSE FIRST DATA POINT')
}


data_list <- build_data(folder = data_folder,
                        dsets = labels$dset_strings_to_use,
                        design_file = design_file,
                        range_file = range_file,
                        
                        #Decisions for every analysis
                        design_names = c('A','B','C','D'),
                        est_err_ell = est_err_ell_run,
                        ell_err = ell_run,
                        
                        #Error handling for Sigma_E
                        split_pb_error = FALSE,
                        full_cor_au = FALSE,
                        
                        #Specifically for MATTER or LBT
                        subset_high_pT_pbpb = subset_for_MATTER,
                        add_header = header_for_LBT,
                        
                        
                        #Things that don't change
                        errors_are_sd = TRUE,
                        write_csvs = FALSE
)


#Specifications of the GPs
R_run <- 5 #Number of PCs
kernel_type_run <- 'matern_5_2' #Nugget used, if not estimating
nugget_est_run <- TRUE #Boolean: estimate nugget? Usually yes
nugget_run = 0 #what covariance kernel

GP_list_cur <- train_GPs(design_input = data_list$design_input,
                         computer_output = data_list$Y,
                         
                         range_file = range_file,
                         
                         
                         #Specifications of the GPs
                         R = R_run, 
                         nugget = nugget_run, 
                         nugget.est = nugget_est_run, 
                         kernel_type = kernel_type_run) 



cal_params_scaled <- calibration_sampler(niters = 1E4,
                                         burnin_prop = 0.3,
                                         t_kap = 1E-3,
                                         
                                         
                                         design_output = data_list$Y,
                                         GP_list = GP_list_cur, #Output from t
                                         y_exp = data_list$y_exp,
                                         Sigma_E = data_list$Sigma_E,
                                         
                                         include_cov_extra = TRUE,
                                         proposal_cor_mat = diag(dim(GP_list_cur[[1]]@input)[2]),
                                         upper_theta_trunc = 1,
                                         bad_corners = list())



## Rescale models ##
params <- rescale_posteriors(cal_params_scaled,
                             ranges = data_list$ranges)

## Heatmaps to check results ##
## See the function for more save options ##
plot_heatmaps(params = params,
              ranges = data_list$ranges,
              pairs = list(c(1,2),c(3,4)),
              save_pics = FALSE,
              design = data_list$design_input,
              plt_points = TRUE)

## Emulator Predictions ##


## Calibration Pairplot ##

#Load all the datasets into a big list
save_path <- 'Results/deliverable_2_6/'

suffixes <- c('m1_all_matern_5_2_5', 'm1_au_matern_5_2_5', 'm1_pb_matern_5_2_5')

param_lists <- vector('list',3)
for(i in 1:3){
  load(paste0(save_path,'save_list_',suffixes[i],'.Rdata'))
  if(i==1){
    ranges_t <- save_list$ranges_t
  }
  param_lists[[i]] <- save_list$transformed_params
}

make_combined_pairplot(param_lists = param_lists, #list of parameter posterior draws
                       labels = c('Combined','RHIC','LHC'), #c('Combined', 'LHC', 'RHIC') usually
                       ranges = ranges_t, #list of parameter ranges
                       save_pics = FALSE,
                       save_path = '', #Default is current directory
                       save_end = '',
                       cols_to_use = 1:ncol(param_lists[[1]]),
                       scatter_alpha = 0.01,
                       num_samples = 1E4, #Number of samples to use.
                       col_names = c('A','C','B','D','Q'),
                       make_pdf = FALSE)


pred_Y <- predict_Y(draws = cal_params_scaled, #matrix of unscaled posterior calibration draws, from calibration sampler
                    GP_list = GP_list_cur, #list of trained GPs
                    design_output = data_list$Y,
                    num_samples = 1E3
)



plot_draws_together(pred_Y$pred_mat,                             
                    
                    dset_labels = labels$systems_to_calibrate,
                    data_orig_labels = labels$data_names_to_use,
                    cen_groups = labels$centrality_groups,
                    cen_labels = labels$centrality_names,
                    
                    include_median_prediction = FALSE,
                    median_prediction_vals = pred_Y$med_pred,
                    
                    original_dset_list = data_list$original_dset_list,
                    
                    
                    col_str_vec = c('red','blue'),
                    col_pts_vec = c('darkred','darkblue'),
                    
                    
                    alpha_val = 0.005,
                    include_legend = TRUE,
                    include_arrows = TRUE,
                    legend_spots = replicate(length(labels$centrality_groups),'topleft'),
                    
                    
                    save_pics = TRUE,
                    title_end = "",
                    save_path = "",
                    save_end = ""
)

