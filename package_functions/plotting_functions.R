## Plotting functions

plot_heatmaps <- function(params,ranges,
                          param_names = colnames(params),
                          pairs = NULL,
                          save_param_names = NULL,
                          save_pics = FALSE,save_path = "", save_end = "", title = 'Posterior Heatmap',
                          plt_points = FALSE,design = NULL, ...){
  if(save_pics & save_path==""){
    print('Warning: saving plots to working directory')
  }
  
  #If you don't specify pairs, just do all of them
  if(is.null(pairs)){
    num_params = dim(params)[2]
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
    f1 <- kde2d(params[,k], (params[,l]), n = 100,
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
          mgp = c(3,1,0), ...
    )
    title(ylab = param_names[l], 
          mgp = c(2.1,1,0),
          cex.lab = 2,
          cex.axis = 1.3)
    
    perc_lvl = c(.6,.75,.9)
    HPDregionplot(params[,c(k,l)], prob = perc_lvl,
                  col=c("black"), lty = c(1,5,3), add=TRUE)
    
    if(plt_points){
      points(design[,k],design[,l],pch = 19, cex = 0.7)
    }

    if(save_pics) dev.off()
  }
}
