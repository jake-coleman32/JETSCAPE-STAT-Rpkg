library(MASS)
library(cluster)
library(emdbook)
library(dplyr)
library(reshape2)
library(RobustGaSP)
library(mvtnorm)
library(Matrix)
library(stringr)
library(GGally)

col_alpha <- function(col_str, alpha){
  if(alpha<0|alpha>1){
    stop('Alpha must be within (0,1)')
  }
  rgb_mat <- col2rgb(col_str,alpha = TRUE)
  return(rgb(rgb_mat[1],rgb_mat[2],rgb_mat[3],alpha = 255*alpha, maxColorValue = 255))
}

# base_save_path = "/Users/Jake/Box Sync/Research/Computer_Emulation/Shanshan_Project/Results/"
# save_folder = ""
# 
# save_path = paste0(base_save_path,save_folder)
##########################################################################
####These are plots to assess the sampler from jetscape_2_concatenate.R
##########################################################################

#Test for Normality
#None are <0.05 even without accounting for multiple testing
normality_test <- function(Y_final){
  for(j in 1:dim(Y_final)[2]){
    print(shapiro.test(Y_final[,j])$p.value)
  }
}

#How many PCs?
calc_Vq <- function(Y_svd,save_pics = FALSE, plotit = FALSE){
  eigs <- Y_svd$d^2
  (V_q <- cumsum(eigs)/sum(eigs))
  
  if(plotit){
    if(save_pics) pdf(paste0(save_path,'var_explained.pdf'))
    plot(V_q[1:6],type = 'o',
         xlab = "Total Number of PCs R",
         main = "Fraction of Variance Explained",
         ylab = expression(F[R]),
         cex.lab = 2,
         cex.axis = 1.3,
         cex.main = 1.4,
         mgp = c(2.3,1,0),
         pch = 19)
    if(save_pics) dev.off()
  }
  return(V_q)
}

###########
##Results
#########
plot_heatmaps <- function(param_plot,ranges,
                          param_names = colnames(param_plot),
                          pairs = NULL,
                          save_param_names = NULL,
                          save_pics = FALSE,save_path = "", save_end = "", title = 'Posterior Heatmap',
                          plt_points = FALSE,design = NULL, ...){
  if(save_pics & save_path==""){
    print('Warning: saving plots to working directory')
  }
  
  #If you don't specify pairs, just do all of them
  if(is.null(pairs)){
    num_params = dim(param_plot)[2]
    num_pairs = sum(1:(num_params-1))
    i = 0
    pairs = vector('list',num_pairs)
    for(k in 1:(num_params-1)){
      for(l in k:num_params)
        i = i + 1
      pairs[[i]] <- c(k,l)
    }
  }
  
  num_pairs = length(pairs)
  
  if(is.null(save_param_names)){
    save_param_names <- vector('list',num_pairs)
  }
  for(i in 1:num_pairs){
    (k = pairs[[i]][1])
    (l = pairs[[i]][2])
    f1 <- kde2d(param_plot[,k], (param_plot[,l]), n = 100,
                lims = c(ranges[[k]],ranges[[l]]))
    
    if(is.null(save_param_names[[i]])){
      save_param_names[[i]] <- paste0(param_names[k],'_',param_names[l])
      print(save_param_names[[i]])
    }
    
    if(save_pics) pdf(paste0(save_path,'heatmap_',save_param_names[[i]],save_end,'.pdf'))
    par(mar=c(5.1, 5.1, 4.1, 2.1))
    image(f1,
          xlab = param_names[k],
          ylab = "",
          main = title,
          cex.lab = 2,
          cex.axis = 1.3,
          cex.main = 1.4,
          mgp = c(3,1,0), ...#,
          # xlim = ranges[[1]],
          # ylim = ranges[[2]]
    )
    title(ylab = param_names[l], 
          mgp = c(2.1,1,0),
          cex.lab = 2,
          cex.axis = 1.3)
    
    perc_lvl = c(.6,.75,.9)
    HPDregionplot(param_plot[,c(k,l)], prob = perc_lvl,
                  col=c("black"), lty = c(1,5,3), add=TRUE)
    
    if(plt_points){
      points(design[,k],design[,l],pch = 19, cex = 0.7)
    }
    # legend('topright',paste0(perc_lvl*100,"%"),title = "Highest Density Kernel Estimate",
    #        lty = c(1,5,3),
    #        bty="n")
    if(save_pics) dev.off()
  }
}

make_param_plot <- function(res_params, ranges){
  param_plot <- matrix(0,dim(res_params)[1],dim(res_params)[2])
  for(j in 1:dim(param_plot)[2]){
    param_plot[,j] <- res_params[,j]*(ranges[[j]][2] - ranges[[j]][1]) + ranges[[j]][1]
  }
  colnames(param_plot) <- names(ranges)
  return(param_plot)
}

unmake_param_plot <- function(param_plot, ranges){
  res_params <- matrix(0,dim(param_plot)[1], dim(param_plot)[2])
  for(j in 1:dim(res_params)[2]){
    res_params[,j] <- (param_plot[,j] - ranges[[j]][1])/(ranges[[j]][2] - ranges[[j]][1])
  }
  colnames(res_params) <- names(ranges)
  return(res_params)
}


transform_draws <- function(param_plot,cols_to_transform = c(1,2),
                            new_names = c('A','C')){
  X = param_plot[,cols_to_transform[1]]
  Y = param_plot[,cols_to_transform[2]]
  new_param_plot = param_plot
  new_param_plot[,1] = X*Y
  new_param_plot[,2] = X*(1-Y)
  colnames(new_param_plot)[cols_to_transform] = new_names
  return(new_param_plot)
}

transform_ranges <- function(ranges,vars_to_transform = c(1,2),
                             new_names = c('A','C')){
  ranges_t <- ranges
  ranges_t[[1]] <- c(0,ranges[[1]][2])
  ranges_t[[2]] <- c(0, ranges[[1]][2])
  names(ranges_t)[1:2] <- new_names
  return(ranges_t)
}

heat_pairplot <- function(param_plot,
                          ranges,
                          split_num = dim(param_plot)[1],
                          param_names = colnames(param_plot), 
                          title = 'Calibration Pairplot',save_pics = FALSE,
                          save_path = "",
                          save_end = "",
                          include_legend = FALSE,
                          legend_names = NULL,
                          legend_cols = NULL,
                          off_diag = "heatmap"){ 
  
  panel.hist <- function(x,col_hist,split_num, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    
    
    h <- hist(x[1:split_num], plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y,xlim = c(0,0.35),col = rgb(1,0,0,.1),  ...)
    
    if(split_num<length(x)){
      h <- hist(x[(split_num+1):length(x)], plot = FALSE)
      breaks <- h$breaks; nB <- length(breaks)
      y <- h$counts; y <- y/max(y)
      rect(breaks[-nB], 0, breaks[-1], y,xlim = c(0,0.35),col = rgb(0,0,1,.1),  ...)
    }
  }
  
  panel.dens <- function(x,split_num,...){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    
    plot(density(x[1:split_num],))
    
    h <- hist(x[1:split_num], plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y,xlim = c(0,0.35),col = rgb(1,0,0,.1),  ...)
    
    if(split_num<length(x)){
      h <- hist(x[(split_num+1):length(x)], plot = FALSE)
      breaks <- h$breaks; nB <- length(breaks)
      y <- h$counts; y <- y/max(y)
      rect(breaks[-nB], 0, breaks[-1], y,xlim = c(0,0.35),col = rgb(0,0,1,.1),  ...)
    }
  }
  
  
  panel.scatter <- function(x,y,split_num,...){
    first_pts = 1:split_num
    
    points(x[first_pts],y[first_pts],#pch = 19,
           col = rgb(1,0,0,.01),cex=.7,
           ...)
    
    #If you actually want to split
    if(split_num<length(x)){
      second_pts <- (split_num+1):length(x)
      points(x[second_pts],y[second_pts],#pch =19,
             col = rgb(0,0,1,.01), cex = 0.7,
             ...)
    }
    
    
  }
  
  
  panel.image <- function(x,y,split_num,...){
    f1 <- kde2d(x, y, n = 100)
    image(f1,add = TRUE)
  }
  
  if(off_diag=="scatter"){
    panel_off = panel.scatter
  }else if(off_diag=="heatmap"){
    if(split_num<dim(param_plot)[1]){
      stop('Heatmap only appropriate for one set of calibration results')
    }
    panel_off = panel.image
  }else{stop("Whoops off_diag parameter is wrong")}
  
  
  
  if(save_pics){ pdf(paste0(save_path,'calibration_pairplot',save_end,".pdf"))}
  pairs(param_plot, 
        panel = panel_off,
        diag.panel = panel.hist,
        #pch = 19,
        #cex = .3,
        # col = rgb(1,0,0,.1),
        labels = param_names,
        upper.panel = NULL,
        cex.lab = 1.5,
        cex.axis = 1.3,
        cex.main = 1.4,
        mgp = c(2.3,1,0),
        #col_hist = 'green',
        las = 2,
        main = title,
        split_num = split_num)
  if(include_legend){
    legend('topright',legend_names,col = legend_cols,
           pch = 19,inset = .13)
  }
  if(save_pics) dev.off()
}



gg_heat_pairs <- function(all_data,
                          ranges,
                          cols_to_use = 1:4,
                          scatter_alpha = 0.1,
                          col_names = c('A','B','C','D')){
  both_col <- 'red3'
  lhc_col <- 'deepskyblue'
  rhic_col <- 'forestgreen'
  
  my_dens <- function(data, mapping, ...) {
    ggplot(data = data, mapping=mapping) +
      geom_density(..., alpha = 0.7, size = 1.1) +
      scale_fill_manual(values=c( both_col, lhc_col, rhic_col )) +
      scale_color_manual(values=c( both_col,lhc_col, rhic_col)) 
  }
  
  
  my_scatter <- function(data,mapping,...){
    ggplot(data = data, mapping = mapping) +
      geom_jitter(alpha = scatter_alpha,size = .3) + 
      guides(col = guide_legend(override.aes = list(shape = 15,
                                                    size = 8,
                                                    alpha = 1),
                                label.theme = element_text(size = 12))) +
      scale_fill_manual(values=c( both_col, lhc_col, rhic_col )) +
      scale_color_manual(values=c( both_col,lhc_col, rhic_col)) + 
      theme(legend.title=element_blank())
  }
  
  p <- ggpairs(all_data, mapping = aes(color = Collider), columns = cols_to_use,
               lower = list(continuous = my_scatter),
               diag = list(continuous = my_dens),
               upper = list(continuous = 'blank'),
               columnLabels = col_names,
               labeller = 'label_parsed'
               #,legend = c(2,1)
  )
  
  
  for(i in 1:length(cols_to_use)){
    for(j in 1:length(cols_to_use)){
      if(j<=i){
        p[i,j] <- p[i,j] + xlim(ranges[[j]][1],ranges[[j]][2]) + 
          theme(axis.text = element_text(size = 11.5))
      }
    }
    
  }
  
  #p <- p + theme(strip.text = element_text(size = 15))
  
  #Remove yaxis on first plot, to remove confusion
  #p[1,1] <- p[1,1] + theme(axis.text.y = element_blank(),
   #                        axis.ticks = element_blank())
  
  
  
  points_legend <- gglegend(my_scatter)
  
  p[1,p$ncol] <- points_legend(all_data, ggplot2::aes(all_data[1,cols_to_use[1]], 
                                                      all_data[1,cols_to_use[2]],
                                                      color = Collider))
  
  #p[1,p$ncol-1] <- points_legend(all_data, ggplot2::aes(all_data[1,cols_to_use[1]], 
   #                                                   all_data[1,cols_to_use[2]],
  #                                                    color = Collider))
  
  
  #print(p + theme(strip.text = element_text(size = 15)))
  print(p)
  
}

make_combined_pairplot <- function(save_lists,
                                   labels,
                                   ranges,
                                   save_pics = FALSE,
                                   save_path = '',
                                   save_end = '',
                                   cols_to_use = 1:5,
                                   scatter_alpha = 0.1,
                                   num_samples = 1E4,
                                   str_param_to_use = 'transformed_params',
                                   col_names = c('A','B','C','D'),
                                   make_pdf = FALSE){
  load(save_lists[1])
  data_to_plot <- matrix(,0,ncol(save_list[[str_param_to_use]]))
  
  if(length(save_lists)!=length(labels)){
    stop('lists and labels not same length')
  }
  for(i in 1:length(save_lists)){
    load(save_lists[i])
    cur_data <- as.data.frame(save_list[[str_param_to_use]]) %>%
      dplyr::sample_n(.,num_samples)
    cur_data$Collider <- labels[i]
    data_to_plot <- rbind(data_to_plot, cur_data)
  }
  
  if(save_pics){
    if(make_pdf){
      pdf(paste0(save_path,'combined_pairplot',save_end,'.pdf'))
    }else{
      png(paste0(save_path,'combined_pairplot',save_end,'.png'))
    }
  }
  gg_heat_pairs(data_to_plot,
                ranges = ranges,
                cols_to_use = cols_to_use,
                col_names = col_names,
                scatter_alpha = scatter_alpha)
  if(save_pics) dev.off()
}


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
                       legend_loc = 'top_left',
                       include_med = FALSE,
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
  if(include_med){
    lines(exp_dat[,pT_col], med_pred_vec,lty = 2, lwd = 2)
  }
  
  
  ##############
  ### DEPRECIATED: NOW WE PUT THE POINTS AND ARROWS IN AFTER ALL LINES
  ###############
  # if(!add){
  #   exp_pch = 15
  # }else{
  #   exp_pch = 17
  # }
  # 
  # ## Add experimental points
  # points(exp_dat[,pT_col],exp_dat[,exp_col],
  #        pch = exp_pch, col = col_pts)
  # 
  # #ylim = c(0,1))
  # if(include_arrows){
  #   arrow_stat_col = col_pts
  #   arrow_sys_col = col_pts
  #   
  #   #Statistical errors
  #   arrows(exp_dat[,pT_col]-arrow_offset, exp_dat[,exp_col] - exp_dat[,err_col_stat],
  #          exp_dat[,pT_col]-arrow_offset, exp_dat[,exp_col] + exp_dat[,err_col_stat],
  #          length=0.05, angle=90, code=3,col = arrow_stat_col)
  #   
  #   #Systematic errors
  #   arrows(exp_dat[,pT_col]+arrow_offset, exp_dat[,exp_col] - exp_dat[,err_col_sys],
  #          exp_dat[,pT_col]+arrow_offset, exp_dat[,exp_col] + exp_dat[,err_col_sys],
  #          length=0.05, angle=90, code=3, col = arrow_sys_col)
  # }
  
  
  
  if(include_legend){
    legend_cex = 1.1
    if(include_med){
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


plot_draws_together <- function(draws,
                                train_mod,
                                rot_mat,
                                original_dset_list,
                                q_index = NA,
                                condition_on_q = FALSE,
                                dset_labels = NULL,
                                data_orig_labels = NULL,
                                cen_groups = NULL,
                                cen_labels = NULL,
                                RAA_med_vals = NULL,
                                med_pred_type = 'RAA',
                                col_str_vec = c('red','blue'),
                                col_pts_vec = c('darkred','darkblue'),
                                include_med = FALSE,
                                scale_mean_vec = NULL, #sd_vec
                                scale_sd_vec = NULL, #y_means
                                delta_list = NULL,
                                num_samples = 1E3,
                                title_end = "",
                                alpha_val = 0.01,
                                save_pics = FALSE,
                                save_path = "",
                                save_end = "",
                                include_legend = TRUE,
                                include_arrows = TRUE,
                                legend_spots = replicate(J,'top_left')){
  
  if(is.null(scale_sd_vec)|is.null(scale_mean_vec)){
    stop('Must set values for scale_sd_vec and scale_mean_vec. Try sd_vec and y_means')
  }else if(dim(draws)[1]<num_samples){
    stop('Num samples bigger than number of draws')
  }else if(sum(do.call(c,lapply(cen_groups,length)) - length(col_str_vec))){
    print('Centrality Groups:')
    print(cen_group)
    print('Color vector:')
    print(col_str_vec)
    stop("Color vec and centrality group lengths not all the same")
  }
  
  Y_pred = exp_predictions(theta = draws,
                           model_list = train_mod,
                           y_sd = scale_sd_vec,
                           y_means = scale_mean_vec,
                           Y_svd = comp_mod$Y_svd,
                           cov_extra = comp_mod$cov_extra,
                           q = length(train_mod))$mean %>% t()
  
  
  Y_pred <- Y_pred[sample(1:dim(Y_pred)[1],num_samples),]
  
  
  if(include_med){
    if(med_pred_type=='emulator'){
      
      if(condition_on_q){
        q_mode <- density(draws[,q_index])$x[which.max(density(draws[,q_index])$y)]
        q_width = 0.3
        good_q <- which(draws[,q_index] > (q_mode - q_width) &
                          draws[,q_index] < (q_mode + q_width) )
        draws <- draws[good_q, ]
      }
      
      (param_meds <- apply(draws,2,median) %>% t())
      
      (med_pred <- exp_predictions(theta = param_meds,
                                   model_list = train_mod,
                                   y_sd = scale_sd_vec,
                                   y_means = scale_mean_vec,
                                   Y_svd = comp_mod$Y_svd,
                                   cov_extra = comp_mod$cov_extra,
                                   q = length(train_mod))$mean %>% t())
    }else if(med_pred_type=='RAA'){
      (med_pred = RAA_med_vals)
    }else{
      stop('med_pred_type must be one of "emulator", "RAA"')
    }
    print(paste0('Median prediction values are ',med_pred_type))
  }else{
    med_pred <- matrix(0,1,length(scale_mean_vec))
  }
  
  cur_col <- 1
  for(j in 1:length(cen_groups)){
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
      cur_dset <- t(Y_pred[,cur_col:(cur_col + num_pT - 1)])
      
      if(med_pred_type=='emulator'){
        (cur_med_pred <- med_pred[,cur_col:(cur_col + num_pT - 1)])
      }else{
        (cur_med_pred <- RAA_med_vals[[j]][[k]])
      }
      
      
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
                   delta = delta_list[[dset_index]],
                   alpha_val = alpha_val,
                   include_legend = include_legend,
                   legend_loc = legend_spots[dset_index],
                   include_med = include_med,
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
                   delta = delta_list[[dset_index]],
                   alpha_val = alpha_val,
                   include_legend = include_legend,
                   legend_loc = legend_spots[dset_index],
                   include_med = include_med,
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


#Predict PCA vals
###DEPRECIATED######
plot_draws <- function(draws,
                       train_mod,
                       rot_mat,
                       original_dset_list,
                       include_med = FALSE,
                       dset_labels = NULL,
                       scale_mean_vec = NULL, #sd_vec
                       scale_sd_vec = NULL, #y_means
                       delta_list = NULL,
                       num_samples = 1E3,
                       title_end = "",
                       alpha_val = 0.01,
                       save_pics = FALSE,
                       save_path = "",
                       save_end = "",
                       include_legend = TRUE,
                       include_arrows = TRUE,
                       legend_spots = replicate(J,'top_left')){
  
  if(is.null(scale_sd_vec)|is.null(scale_mean_vec)){
    stop('Must set values for scale_sd_vec and scale_mean_vec. Try sd_vec and y_means')
  }else if(dim(draws)[1]<num_samples){
    stop('Num samples bigger than number of draws')
  }
  
  Y_pred = exp_predictions(theta = draws,
                           model_list = train_mod,
                           y_sd = scale_sd_vec,
                           y_means = scale_mean_vec,
                           Y_svd = comp_mod$Y_svd,
                           cov_extra = comp_mod$cov_extra,
                           q = length(train_mod))$mean %>% t() #%>%
  #dplyr::slice(sample(1:dim(Y_pred)[1],num_samples))
  
  
  
  # ##Predict the mu_Z, collect them
  # colnames(draws) <- colnames(ranges)
  # post_pred_mod <-lapply(train_mod, RobustGaSP::predict,
  #                        #  testing_trend = as.data.frame(cbind(rep(1,dim(draws)[1]),draws)),
  #                        testing_input = as.data.frame(draws))
  # 
  # pred_means <- lapply(post_pred_mod,function(x)x$mean) %>%
  #   do.call(cbind,.)
  # 
  # 
  # #Rotate
  # Y_pred_scaled <- pred_means %*%t(rot_mat)
  # 
  # #Scale
  # 
  # Y_pred <- sweep(Y_pred_scaled,2,scale_sd_vec,FUN = "*") %>%
  #   sweep(2,scale_mean_vec,FUN = '+') 
  # #Y_pred <- sweep(Y_pred_scaled,2,apply(scale_mat,2,mean),FUN = '+') 
  # 
  
  
  Y_pred <- Y_pred[sample(1:dim(Y_pred)[1],num_samples),]
  
  
  if(include_med){
    (param_meds <- apply(draws,2,median) %>% t())
    
    (med_pred <- exp_predictions(theta = param_meds,
                                 model_list = train_mod,
                                 y_sd = scale_sd_vec,
                                 y_means = scale_mean_vec,
                                 Y_svd = comp_mod$Y_svd,
                                 cov_extra = comp_mod$cov_extra,
                                 q = length(train_mod))$mean %>% t())
    
  }else{
    med_pred <- matrix(0,1,length(scale_mean_vec))
  }
  
  cur_col <- 1
  for(j in 1:length(original_dset_list)){
    (num_pT <- dim(original_dset_list[[j]])[1])
    cur_dset <- t(Y_pred[,cur_col:(cur_col + num_pT - 1)])
    
    cur_med_pred <- med_pred[,cur_col:(cur_col + num_pT - 1)]
    
    exp_dset <- cbind('pT' = original_dset_list[[j]]$pT,
                      'RAA_exp' = original_dset_list[[j]]$RAA_exp,
                      'exp_err_stat' = original_dset_list[[j]]$Stat_err,
                      'exp_err_sys' =  original_dset_list[[j]]$Sys_err)
    (cur_col = cur_col + num_pT)
    
    
    #Plot 'em
    if(save_pics) pdf(paste0(save_path,dset_labels[j],save_end,'.pdf'))
    plot_lines(pred_dat = cur_dset,
               exp_dat = exp_dset,
               med_pred_vec = cur_med_pred,
               delta = delta_list[[j]],
               main = paste(dset_labels[j],title_end),
               alpha_val = alpha_val,
               include_legend = include_legend,
               legend_loc = legend_spots[j],
               include_med = include_med)
    if(save_pics) dev.off()
  }
  
}


plot_discrep_fun <- function(discrep_list,
                             ell_mat,
                             lambda_mat,
                             pT_train,
                             pT_pred,
                             alpha = 2,
                             include_pts = TRUE,
                             dset_labels = dset_strings,
                             save_pics = FALSE,
                             save_path = "",
                             save_suffix = ""){
  
  mu_delt <- sig_delt <- vector('list',length(discrep_list))
  for(j in 1:length(discrep_list)){
    y = discrep_list[[j]]
    x_train = matrix(pT_train[[j]])
    x_pred = matrix(sort(c(pT_train[[j]],pT_pred[[j]])))
    ell_j = mean(ell_mat[,j])
    lambda_j = mean(lambda_mat[,j])
    sig_22 = cov_mat_multi(x_train,
                           x_train,
                           ell_vec = ell_j,
                           lambda = lambda_j,
                           alpha_vec = alpha,
                           nugget = 1E-6)
    sig_22_inv <- chol2inv(chol(sig_22))
    sig_12 = cov_mat_multi(x_pred,
                           x_train,
                           ell_vec = ell_j,
                           lambda = lambda_j,
                           alpha_vec = alpha,
                           nugget = 1E-6)
    sig_21 = t(sig_12)
    sig_11 = cov_mat_multi(x_pred,
                           x_pred,
                           ell_vec = ell_j,
                           lambda = lambda_j,
                           alpha_vec = alpha,
                           nugget = 1E-6)
    
    mu_delt[[j]] <- sig_12%*%sig_22_inv%*%y
    sig_delt[[j]] <- diag(sig_11 - sig_12%*%sig_22_inv%*%sig_21)
    if(save_pics)pdf(paste0(save_path,"discrep_",dset_labels[j],save_suffix,".pdf"))
    plot(x_pred,mu_delt[[j]],type = 'l',
         # ylim = c(min(mu_delt[[j]] - 2*sqrt(sig_delt[[j]])),
         #          max(mu_delt[[j]] + 2*sqrt(sig_delt[[j]]))))
         ylim = c(-.23,.23),
         xlab = expression(p[T]~"(GeV)"),
         #ylab = expression(delta[j]),
         main = paste0("Discrepancy for ",dset_labels[j]),
         cex.lab = 2,
         cex.main = 2,
         ylab = "")
    title(ylab = expression(delta[j]),cex.lab = 2,mgp=c(2,1,0))
    lines(x_pred,mu_delt[[j]] - 2*sqrt(sig_delt[[j]]),col = 'blue',
          lty = 2)
    lines(x_pred,mu_delt[[j]] + 2*sqrt(sig_delt[[j]]),col = 'blue',
          lty = 2)
    abline(h=0)
    if(include_pts){
      points(x_train,y)
    }
    if(save_pics) dev.off()
  }
  invisible(list(mu_delt = mu_delt,sig_delt = sig_delt))
}

###########
##Validation Plot
##########

##Requires:
#Holdout point
#Design 
#Train dset
# V
# y_sd
# y_mean

##Ugh...should change this to use the new prediction function
predict_holdout <- function(holdout,
                            design = comp_mod$all_data$scaled_d,
                            Y_scaled = comp_mod$Y_final, 
                            #Y_orig = Y, 
                            y_sd_cur = comp_mod$y_sd,
                            y_mean_cur = comp_mod$y_means,
                            q = length(comp_mod$train_mods), plotit = FALSE,
                            verbose = FALSE,
                            cover_level = 95){
  train_d <- design[-holdout, ]
  test_d <- design[holdout, ]
  
  train_Y <- Y_scaled[-holdout, ]
  
  Y_orig <- sweep(Y_scaled,2,y_sd_cur,"*")%>%
    sweep(2,y_mean_cur, "+")
  
  test_Y <- Y_orig[holdout, ]
  
  V_train <- svd(train_Y)$v[ , 1:q]
  
  train_Z = as.matrix(train_Y)%*%V_train
  
  sink("NUL")
  train_mod <- lapply(as.data.frame(train_Z),rgasp,design = train_d,
                      kernel_type = 'pow_exp')
  sink()
  
  test_mod <- lapply(train_mod,RobustGaSP::predict,testing_input = test_d)
  
  pred_Z <- lapply(test_mod,function(x)x$mean) %>%
    do.call(c,.) %>%
    t()
  
  pred_err_Z <- lapply(test_mod,function(x)x$sd^2) %>%
    do.call(c,.) 
  
  pred_Y <- pred_Z %*% t(V_train) %*% diag(y_sd_cur) + y_mean_cur
  
  pred_err_Y <- diag(y_sd_cur) %*% V_train %*% diag(pred_err_Z) %*% t(V_train) %*% diag(y_sd_cur) %>%
    diag()
  
  
  if(cover_level == 95){
    sd_mult = 2
  }else if(cover_level == 60){
    sd_mult = 1
  }else{
    stop('Must pick 95% or 60% for cover level')
  }
  
  
  if(plotit){
    plot(as.numeric(test_Y),as.numeric(pred_Y),pch = 19,cex = .5,
         xlab = 'Holdout Values',
         ylab = 'Predicted Values',
         main = paste0('Emulator Prediction, Holdout ',holdout),
         cex.lab = 2,
         cex.axis = 1.3,
         cex.main = 1.4,
         mgp = c(2.4,1,0),
         ylim = c(min(as.numeric(pred_Y) - sd_mult*sqrt(pred_err_Y)),
                  max(as.numeric(pred_Y) + sd_mult*sqrt(pred_err_Y))))
    arrows(as.numeric(test_Y), as.numeric(pred_Y) - sd_mult*sqrt(pred_err_Y),
           as.numeric(test_Y), as.numeric(pred_Y) + sd_mult*sqrt(pred_err_Y),
           length=0, angle=90, code=3)
    
    abline(a = 0,b = 1)
  }
  
  good_pred_95 <- test_Y < as.numeric(pred_Y) + 2*sqrt(pred_err_Y) & 
    test_Y > as.numeric(pred_Y) - 2*sqrt(pred_err_Y)
  
  good_pred_60 <- test_Y < as.numeric(pred_Y) + sqrt(pred_err_Y) & 
    test_Y > as.numeric(pred_Y) - sqrt(pred_err_Y)
  
  
  return(list(mean_cover_95 = mean(good_pred_95), mean_cover_60 = mean(good_pred_60),
              mean_width = mean(4*sqrt(pred_err_Y))))
}


collect_all_holdouts <- function(D_c_theta = comp_mod$all_data$scaled_d,
                                 Y_scaled = comp_mod$Y_final, 
                                 #Y_orig = Y, 
                                 y_sd_cur = comp_mod$y_sd,
                                 y_mean_cur = comp_mod$y_means,
                                 q = length(comp_mod$train_mods)){
  
  n_holdouts <- nrow(D_c_theta)
  cover_95 <- cover_60 <- widths <-  numeric(n_holdouts)
  for(h in 1:n_holdouts){
    invisible(pred_hold_h <- predict_holdout(holdout = h,
                                             design = D_c_theta,
                                             Y_scaled = Y_scaled,
                                             y_sd_cur = y_sd_cur,
                                             y_mean_cur = y_mean_cur,
                                             q = q))
    cover_95[h] <- pred_hold_h$mean_cover_95
    cover_60[h] <- pred_hold_h$mean_cover_60
    widths[h] <- pred_hold_h$mean_width
    
    flush.console()
    cat("\r FINISHED HOLDOUT ", h ,'/',n_holdouts)
    
  }
  return(list(cover_95 = cover_95, cover_60 = cover_60, widths = widths))
}


