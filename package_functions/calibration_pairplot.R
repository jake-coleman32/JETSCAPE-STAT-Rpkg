gg_heat_pairs <- function(all_data, #Data frame of posterior draws
                          ranges, #List of parameter ranges
                          cols_to_use = 1:ncol(all_data), #In case you want to choose specific parameters
                          scatter_alpha = 0.1, # controls the shade of the scatter plot points. Larger is darker
                          col_names = c('A','B','C','D')){
  both_col <- 'red3'
  lhc_col <- 'deepskyblue'
  rhic_col <- 'forestgreen'
  
  #Density function for the diagonal elements
  my_dens <- function(data, mapping, ...) {
    ggplot(data = data, mapping=mapping) +
      geom_density(..., alpha = 0.7, size = 1.1) +
      scale_fill_manual(values=c( both_col, lhc_col, rhic_col )) +
      scale_color_manual(values=c( both_col,lhc_col, rhic_col)) 
  }
  
  #Scatter plot function for the off-diagonal elements
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
               lower = list(continuous = my_scatter), #lower off-diagonal is scatter plot
               diag = list(continuous = my_dens), #diagonal is density
               upper = list(continuous = 'blank'),#upper off-diagonal is blank
               columnLabels = col_names,
               labeller = 'label_parsed'
  )
  
  #Set the ranges of the x-axis, so that it shows the distribution over the range of the input
  for(i in 1:length(cols_to_use)){
    for(j in 1:length(cols_to_use)){
      if(j<=i){
        p[i,j] <- p[i,j] + xlim(ranges[[j]][1],ranges[[j]][2]) + 
          theme(axis.text = element_text(size = 11.5))
      }
    }
    
  }
  

  #Remove yaxis on first plot, to remove confusion
  p[1,1] <- p[1,1] + theme(axis.text.y = element_blank(),
                         axis.ticks = element_blank())
  
  
  
  points_legend <- gglegend(my_scatter)
  
  p[1,p$ncol] <- points_legend(all_data, ggplot2::aes(all_data[1,cols_to_use[1]], 
                                                      all_data[1,cols_to_use[2]],
                                                      color = Collider))
  
  
  ## This doesn't work for some reason, if you change the text size you remove the legend
  ## Has to do with themes conflicting, I think
  #print(p + theme(strip.text = element_text(size = 15)))
  
  print(p)
  
}

make_combined_pairplot <- function(param_lists, #list of parameter posterior draws
                                   labels, #c('Combined', 'LHC', 'RHIC') usually
                                   ranges, #list of parameter ranges
                                   save_pics = FALSE,
                                   save_path = '', #Default is current directory
                                   save_end = '', #If doing multiple runs or something
                                   cols_to_use = 1:ncol(param_lists[[1]]),
                                   scatter_alpha = 0.1,
                                   num_samples = 1E4, #Number of samples to use.
                                   col_names = c('A','B','C','D'),
                                   make_pdf = FALSE){
  
  #Initially data matrix
  data_to_plot <- matrix(,0,ncol(param_lists[[1]]))
  
  #Sense check
  if(length(param_lists)!=length(labels)){
    stop('lists and labels not same length')
  }
  
  #Concatenating all the different posterior draws
  #Labeling each one with the specified label
  for(i in 1:length(param_lists)){
    cur_data <- as.data.frame(param_lists[[i]]) %>%
      dplyr::sample_n(.,num_samples)
    cur_data$Collider <- labels[i]
    data_to_plot <- rbind(data_to_plot, cur_data)
  }
  
  if(save_pics){
    if(make_pdf){#Makes a very large file
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