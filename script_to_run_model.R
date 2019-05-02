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

data_folder = "Data/to_Jake/MATTER/forSTAT-MATTER/"

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
(dsets_to_use = which(str_detect(all_dsets, paste(systems_to_calibrate,collapse = "|"))))

cur_dsets = all_dsets[dsets_to_use]

system(paste0('ls ',data_folder))
design_file = paste0(data_folder,'matter_lhc.txt')
range_file = paste0(data_folder,'ranges.Rdata')

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


data_list <- build_data(folder = data_folder,
                        dsets = cur_dsets,
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

design_input_matrix <- 
  
  #Specifications of the GPs
  R_run <- 3 #Number of PCs
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



cal_params <- calibration_sampler(niters = 1E4,
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
