## Emulator Prediction Plots

#Predict model output given calibration draws, trained GPs, and original design output
predict_Y <- function(draws, #matrix of unscaled posterior calibration draws, from calibration sampler
                      GP_list, #list of trained GPs
                      design_output, #original computer model output for design
                      num_samples = 1E4){
  
  if(num_samples > dim(draws)[1]){
    stop('Cannot have more samples than posterior draws')
  }
  
  #Renaming for ease of use
  Y <- design_output
  
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  
  y_sd <- apply(Y, 2, sd)
  y_means <- apply(Y, 2, mean)
  
  Y_final <- as.matrix(sweep(Y, 2, y_means, FUN = '-')) %*% diag(1/y_sd)
  Y_svd <- svd(Y_final)
  
  R <- length(GP_list)
  
  
  V_b = Y_svd$v[ ,(R+1):q]
  S_b = diag(Y_svd$d[(R+1):q])
  cov_extra <- 1/(n-1) * diag(y_sd) %*% V_b %*% S_b^2 %*% t(V_b) %*% diag(y_sd)
  
  
  Y_pred = exp_predictions(theta = draws,
                           GP_list = GP_list,
                           y_sd = y_sd,
                           y_means = y_means,
                           Y_svd = Y_svd,
                           cov_extra = cov_extra)$mean %>% t()
  
  
  Y_pred <- Y_pred[sample(1:dim(Y_pred)[1],num_samples),]  
  
  (param_meds <- apply(draws,2,median) %>% t())
  
  (med_pred <- exp_predictions(theta = param_meds,
                               GP_list = GP_list,
                               y_sd = y_sd,
                               y_means = y_means,
                               Y_svd = Y_svd,
                               cov_extra = cov_extra)$mean %>% t())
  
  return(list(pred_mat = Y_pred,
              med_pred = med_pred))
  
}


#Making each individual plot
plot_lines <- function(pred_dat,exp_dat,
                       med_pred_vec = NULL,
                       pT_col = 'pT',
                       exp_col = 'RAA_exp',
                       err_col_stat = 'exp_err_stat', 
                       err_col_sys = 'exp_err_sys',
                       col_str = 'red',
                       col_pts = 'darkred',
                       all_col_pts = c('darkred','darkblue'),
                       all_col_strs = c('red','blue'),
                       alpha_val = 0.01,
                       delta = rep(0,dim(exp_dat)[1]),
                       include_legend = TRUE,
                       legend_loc = 'topleft',
                       include_median_prediction = FALSE,
                       add = FALSE,
                       include_arrows = TRUE,
                       cen_labels = NULL,
                       plot_title = "",
                       data_loc_name = "",
                       ...){
  
  if(!add){
    #xlimits
    pt_lim = min(exp_dat[,pT_col])
    if(str_detect(plot_title,'AuAu200')){
      pt_lim = c(pt_lim,20)
      txt_x = 11.7
    }else if(str_detect(plot_title,'PbPb2760')){
      pt_lim = c(pt_lim,90)
      #txt_x = 45 #for MATTER (?why)
      txt_x = 27
    }else if(str_detect(plot_title,'PbPb5020')){
      pt_lim = c(pt_lim,100)
      #txt_x = 45 #for MATTER (?why)
      txt_x = 30
    }else{
      print('Warning: Using auto-generated tick marks for pT')
      pt_lim = c(pt_lim,max(exp_dat[,pT_col]))
    }
    
    plot(exp_dat[,pT_col],pred_dat[,1],
         cex = 0.5,
         col = col_alpha(col_str,alpha_val),
         type = 'l',
         # ylim = c(min(exp_dat[,exp_col] - 2*exp_dat[,err_col_sys],
         #              exp_dat[,exp_col] - 2*exp_dat[,err_col_stat]),
         #          max(exp_dat[,exp_col] + 2*exp_dat[,err_col_sys],
         #              exp_dat[,exp_col] + 2*exp_dat[,err_col_stat])),
         ylim = c(0,1.8),
         xlim = pt_lim,
         xlab = expression(p[T]~"(GeV)"),
         ylab = expression(R[AA]),
         cex.lab = 2,
         mgp = c(2.3,1,0),
         pch = 15,
         yaxt = 'n',
         xaxt = 'n',
         ...)
    title(main = plot_title, line = -1.5, cex.main = 1.3, adj = .97)
    box(lwd = 3)
    
    #Add data location labels
    text(txt_x, 1.47, paste0('Data from ',data_loc_name),cex = 1, font = 2)
    
    #y axis ticks and labels
    axis(side = 2, at = seq(0,1.5,by = 0.3), las = 1, lwd.ticks = 3, labels = NA)
    axis(side = 2, at = seq(0,1.5,by = 0.3), las = 1, 
         cex.axis = 1.2, lwd = 0, line = -0.4)
    
    #x axis ticks and labels
    if(str_detect(plot_title,'AuAu')){
      pt_ticks = seq(10,20,by = 2)
    }else if(str_detect(plot_title,'PbPb')){
      pt_ticks = seq(10,100, by = 10)
    }else{
      print('Warning: Using auto-generated tick marks for pT')
      pt_ticks = seq(min(exp_dat[,pT_col]),max(exp_dat[,pT_col]), length.out = 5)
    }
    
    axis(side = 1, at = pt_ticks, lwd.ticks = 3, labels = NA)
    axis(side = 1, at = pt_ticks, cex.axis = 1.2, lwd = 0, line = -0.4)
    
    
  }else{
    lines(exp_dat[,pT_col],pred_dat[,1],
          cex = 0.5,
          col = col_alpha(col_str,alpha_val),type = 'l')
  }
  #x_range = max(exp_dat[,pT_col]) - min(exp_dat[,pT_col])
  #arrow_offset = x_range/130
  arrow_offset = 0
  
  for(i in 2:dim(pred_dat)[2]){
    lines(exp_dat[,pT_col],pred_dat[,i] + delta,cex = 0.5,
          col = col_alpha(col_str,alpha_val))
  }
  if(include_median_prediction){
    lines(exp_dat[,pT_col], med_pred_vec,lty = 2, lwd = 2)
  }
  
  
  
  if(include_legend){
    legend_cex = 1.1
    if(include_median_prediction){
      legend(legend_loc, c(paste(cen_labels,'Centrality'),'Median Predictions'),
             lty = c(1,1,2),lwd = c(2,2,2),col = c(all_col_strs,'black'),
             box.lwd = 2, cex = legend_cex)
      
      #Put the point shapes on top
      legend(legend_loc, c(rep("",length(cen_labels)),""), col = all_col_pts,
             pch = c(15,17,NA), bty='n', cex = legend_cex)
      
    }else{
      legend(legend_loc, paste(cen_labels,'Centrality'),lwd = 2, col = all_col_strs,
             box.lwd = 2, cex = legend_cex)
      
      #Put the point shapes on top
      legend(legend_loc, rep("",length(cen_labels)), col = all_col_pts,
             pch = c(15, 17), bty='n', box.lwd = 2, cex = legend_cex)
      
    }
  }
  
}




#Making multiple plots with multiple centralities
plot_draws_together <- function(pred_Y,                             
                                
                                dset_labels = NULL,
                                data_orig_labels = NULL,
                                cen_groups = NULL,
                                cen_labels = NULL,
                                
                                include_median_prediction = FALSE,
                                median_prediction_vals = NULL,
                                
                                original_dset_list = NULL,
                                
                                
                                col_str_vec = c('red','blue'),
                                col_pts_vec = c('darkred','darkblue'),
                                
                                
                                alpha_val = 0.01,
                                include_legend = TRUE,
                                include_arrows = TRUE,
                                legend_spots = replicate(length(cen_groups),'topleft'),
                                
                                
                                save_pics = FALSE,
                                title_end = "",
                                save_path = "",
                                save_end = ""
){
  
  if(sum(do.call(c,lapply(cen_groups,length)) - length(col_str_vec))){
    print('Centrality Groups:')
    print(cen_group)
    print('Color vector:')
    print(col_str_vec)
    stop("Color vec and centrality group lengths not all the same")
  }
  
  if(is.null(original_dset_list)){
    stop('Must include original_dset_list')
  }
  
  
  if(include_median_prediction){
    if(is.null(median_prediction_vals)){
      stop('Cannot include median prediction values if input is NULL')
    }
  }else{
    median_prediction_vals <- matrix(0,1,dim(pred_Y)[2])
  }
  
  cur_col <- 1
  for(j in 1:length(cen_groups)){
    print(j)
    (cur_cen_labels = cen_labels[[j]])
    (cur_group = names(cen_groups)[j])
    (K_j <- length(cen_groups[[j]]))
    (cur_data_loc_name = data_orig_labels[j])
    
    exp_dset_list_j = vector('list',K_j)
    
    for(k in 1:K_j){
      (dset_index = cen_groups[[j]][k])
      (cur_color <- col_str_vec[k])
      (cur_pts_col <- col_pts_vec[k])
      
      
      (num_pT <- dim(original_dset_list[[dset_index]])[1])
      cur_dset <- t(pred_Y[,cur_col:(cur_col + num_pT - 1)])
      
      
      (cur_med_pred <- median_prediction_vals[,cur_col:(cur_col + num_pT - 1)])
      
      
      
      exp_dset <- cbind('pT' = original_dset_list[[dset_index]]$pT,
                        'RAA_exp' = original_dset_list[[dset_index]]$RAA_exp,
                        'exp_err_stat' = original_dset_list[[dset_index]]$Stat_err,
                        'exp_err_sys' =  original_dset_list[[dset_index]]$Sys_err)
      
      exp_dset_list_j[[k]] <- exp_dset
      (cur_col = cur_col + num_pT)
      #Plot 'em
      if(k==1){
        if(save_pics) pdf(paste0(save_path,dset_labels[j],save_end,'.pdf'))
        
        plot_lines(pred_dat = cur_dset,
                   exp_dat = exp_dset,
                   med_pred_vec = cur_med_pred,
                   alpha_val = alpha_val,
                   include_legend = include_legend,
                   legend_loc = legend_spots[j],
                   include_median_prediction = include_median_prediction,
                   col_str = cur_color,
                   all_col_strs = col_str_vec,
                   col_pts = cur_pts_col,
                   all_col_pts = col_pts_vec,
                   cen_labels = cur_cen_labels,
                   include_arrows = include_arrows,
                   
                   plot_title = paste(cur_group,title_end),
                   data_loc_name = cur_data_loc_name
        )
      }else{
        plot_lines(pred_dat = cur_dset,
                   exp_dat = exp_dset,
                   med_pred_vec = cur_med_pred,
                   alpha_val = alpha_val,
                   include_legend = include_legend,
                   legend_loc = legend_spots[j],
                   include_median_prediction = include_median_prediction,
                   col_str = cur_color,
                   all_col_strs = col_str_vec,
                   col_pts = cur_pts_col,
                   all_col_pts = col_pts_vec,
                   cen_labels = cur_cen_labels,
                   include_arrows = include_arrows,
                   add = TRUE)
      }
      if(k==K_j){
        #########
        ## Plot all the points on top of the lines
        #########
        pT_col = 'pT'
        exp_col = 'RAA_exp'
        err_col_stat = 'exp_err_stat'
        err_col_sys = 'exp_err_sys'
        
        #Loop over the datasets in this collision system
        for(i in 1:K_j){
          (col_pts <- col_pts_vec[i])
          
          exp_dat = exp_dset_list_j[[i]]
          if(i < K_j){
            exp_pch = 15
          }else{
            exp_pch = 17
          }
          ## Add experimental points
          points(exp_dat[,pT_col],exp_dat[,exp_col],
                 pch = exp_pch, col = col_pts)
          
          #ylim = c(0,1))
          
          if(include_arrows){
            arrow_stat_col = col_pts
            arrow_sys_col = col_pts
            
            #x_range = max(exp_dat[,pT_col]) - min(exp_dat[,pT_col])
            #arrow_offset = x_range/130
            arrow_offset = 0
            
            #Statistical errors
            arrows(exp_dat[,pT_col]-arrow_offset, exp_dat[,exp_col] - exp_dat[,err_col_stat],
                   exp_dat[,pT_col]-arrow_offset, exp_dat[,exp_col] + exp_dat[,err_col_stat],
                   length=0.05, angle=90, code=3,col = arrow_stat_col)
            
            #Systematic errors
            arrows(exp_dat[,pT_col]+arrow_offset, exp_dat[,exp_col] - exp_dat[,err_col_sys],
                   exp_dat[,pT_col]+arrow_offset, exp_dat[,exp_col] + exp_dat[,err_col_sys],
                   length=0.05, angle=90, code=3, col = arrow_sys_col)
          }
        }
        
        if(save_pics) dev.off()
      }
    }#k loop
  }#j loop
  
}
