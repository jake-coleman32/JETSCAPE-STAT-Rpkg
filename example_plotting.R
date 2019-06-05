library(dplyr)
library(GGally)

package_dir <- '/Users/Jake/Box/Research/Computer_Emulation/Shanshan_Project/JETSCAPE_STAT_pkg/'
setwd(package_dir)

source('package_functions/helper_functions.R')
source('package_functions/calibration_pairplot.R')
source('package_functions/plot_predictions.R')
source('package_functions/exp_predictions.R')


###############
## Calibration Pairplot 
################

## Load LBT+MATTER_1 calibration posteriors, for Combined, LBT, and RHIC
param_lists <- vector('list',3)

param_lists[[1]] <- read.csv('MATTER+LBT_1/posterior_draws_m1_all_matern_5_2_5.csv')
param_lists[[2]] <- read.csv('MATTER+LBT_1/posterior_draws_m1_au_matern_5_2_5.csv')
param_lists[[3]] <- read.csv('MATTER+LBT_1/posterior_draws_m1_pb_matern_5_2_5.csv')

param_ranges <- list('A' = c(0, 1.5),
                     'C' = c(0, 1.5),
                     'B' = c(0, 20),
                     'D' = c(0, 20),
                     'Q' = c(1, 4))


make_combined_pairplot(param_lists = param_lists,
                       labels = c('Combined','RHIC','LHC'), 
                       ranges = param_ranges, 
                       save_pics = TRUE,
                       save_path = 'MATTER+LBT_1/', 
                       save_end = '', 
                       cols_to_use = 1:ncol(param_lists[[1]]),
                       scatter_alpha = 0.01,
                       num_samples = 1E4, 
                       col_names = c('A','C','B','D','Q'),
                       make_pdf = FALSE)


#################
## Prediction Plots
##################

pred_emulator <- read.csv('MATTER+LBT_1/posterior_emulator_predictions.csv')


#Labels to make
# collision system labels
collision_systems <- c('AuAu200','PbPb2760','PbPb5020')
# data origin labels 
data_origins <- c('PHENIX', 'ATLAS', 'CMS')
# centrality groups
centrality_groups <- list('AuAu200' = c(1,2),
                          'PbPb2760' = c(3,4),
                          'PbPb5020' = c(5,6))
# centrality labels
centrality_labels = list('AuAu200' = c("0% - 10%", "40% - 50%"),
                         'PbPb2760' = c("0% - 5%", "30% - 40%"),
                         'PbPb5020' = c("0% - 10%", "30% - 50%"))

#Other things needed
#original data list
load('MATTER+LBT_1/original_dset_list.Rdata')

med_vals_vec <- read.csv('MATTER+LBT_1/median_prediction_values.csv')

plot_draws_together(pred_Y = pred_emulator,                             
                    
                    dset_labels = collision_systems,
                    data_orig_labels = data_origins,
                    cen_groups = centrality_groups,
                    cen_labels = centrality_labels,
                    
                    include_median_prediction = TRUE,
                    median_prediction_vals = med_vals_vec,
                    
                    original_dset_list = original_dset_list,
                    
                    
                    col_str_vec = c('skyblue','pink'),
                    col_pts_vec = c('darkblue','darkred'),
                    
                    
                    alpha_val = 0.05,
                    include_legend = TRUE,
                    include_arrows = TRUE,
                    legend_spots = replicate(length(labels$centrality_groups),'topleft'),
                    
                    
                    save_pics = TRUE,
                    title_end = "",
                    save_path = "MATTER+LBT_1/",
                    save_end = "_posterior_predictions"
)
