##
#Put new data and old data together
##what does this entail?
##Maybe make function that takes most recent all_data folder, new folder, and combines

make_global <- function(lst){
  list2env(lst,globalenv())
}

setwd('/Users/Jake/Box/Research/Computer_Emulation/Shanshan_Project/Data/')
#setwd('Documents/Shanshan_Project/Data/')

View(read.table('extended_design/forSTAT-qhat-4D-add2/PbPb2760-cen-00-05.dat',header = FALSE))

combine_data <- function(new_folder,old_folder,
                         datasets,
                         new_lhc, old_lhc,
                         new_ranges_file, old_ranges_file,
                         write_folder,
                         old_header = TRUE,
                         new_header = FALSE,
                         lhc_name = 'combined_lhc.dat',
                         scaled_d_name = 'combined_scaled_d.dat',
                          write = FALSE){
  for(d in 1:length(datasets)){
    (dset = datasets[d])
    old_dset = read.table(paste0(old_folder, dset),header = old_header)
    new_dset = read.table(paste0(new_folder, dset),header = new_header)
    
    combined_dset = cbind(old_dset,new_dset[,-c(1:4)])
    
    
    if(write){
      write.table(combined_dset,file = paste0(write_folder,dset),
                  col.names = FALSE, row.names = FALSE)
    }
  }
  
  #Annoying: I save the ranges as .Rdata files and they're always called 'ranges'
  load(paste0(new_folder,new_ranges_file))
  new_ranges = ranges
  load(paste0(old_folder,old_ranges_file))
  old_ranges = ranges
  
  #Take the old lhc scaled to actual numbers (i.e. that I gave Shanshan), concatenate
  #Recale by combined ranges
  
  #Make the ranges the 
  combined_ranges = vector('list',length(old_ranges))
  for(i in 1:length(combined_ranges)){
    combined_ranges[[i]] = c(min(old_ranges[[i]][1],new_ranges[[i]][1]),
                             max(old_ranges[[i]][2],new_ranges[[i]][2]))
  }
  
  old_design = read.table(paste0(old_folder, old_lhc),header = TRUE)
  new_design = read.table(paste0(new_folder, new_lhc),header = TRUE)
  
  combined_design = rbind(old_design,new_design)
  
  combined_scaled_d = matrix(0,dim(combined_design)[1],dim(combined_design)[2])
  for(i in 1:dim(combined_design)[2]){
    combined_scaled_d[,i] = (combined_design[,i] - combined_ranges[[i]][1])/(combined_ranges[[i]][2] - combined_ranges[[i]][1])
  }
  
  if(write){
    write.table(combined_design,file = paste0(write_folder,lhc_name))
    write.table(combined_scaled_d,file = paste0(write_folder,scaled_d_name))
    ranges <- combined_ranges
    save(ranges,file = paste0(write_folder,'ranges.Rdata'))
  }
  
}

all_dsets <- c(
  "AuAu200-cen-00-10"
  ,"AuAu200-cen-40-50"
  ,"PbPb2760-cen-00-05"
  ,"PbPb2760-cen-30-40"
  ,"PbPb5020-cen-00-10"
  ,"PbPb5020-cen-30-50"
)

####NOTE: LHC's AND RANGES ARE WRONG HERE BC THE NEW FOLDER DOESN'T HAVE THEM FOR SOME REASON
combine_data(new_folder = 'extended_design/forSTAT-qhat-4D-add2/',
             old_folder = 'discrepancy/forSTAT-combine/',
             datasets = paste0(all_dsets,'.dat'),
             new_lhc = 'discrepancy/latin_hc_design_add_2.txt',
             old_lhc = 'discrepancy/forSTAT-combine/latin_hc_design_added.txt',
             new_ranges_file = 'discrepancy/add_2_ranges.Rdata',
             old_ranges_file = 'discrepancy/forSTAT-combine/add_ranges.Rdata',
             write_folder = 'extended_design/only_shortened/',
             old_header = FALSE,
             new_header = FALSE,
             write = TRUE)

View(read.table('extended_design/all_combined/combined_lhc.dat',header = TRUE))
load('extended_design/all_combined/ranges.Rdata')


combine_data(new_folder =  'five_params/forSTAT-5D-10_more_points/',
             old_folder = 'five_params/forSTAT-5D-const-Q0/',
             datasets = paste0(all_dsets,'.dat'),
             new_lhc = 'latin_hc_five_add.txt',
             old_lhc = 'latin_hc_five.txt',
             new_ranges_file = 'add_ranges.Rdata',
             old_ranges_file = 'ranges.Rdata',
             write_folder = 'five_params/forSTAT-5D-const-Q0-extended/',
             lhc_name = 'const_combined_lhc.txt',
             scaled_d_name = 'const_combined_scaled_d.txt',
             old_header = FALSE,
             new_header = FALSE,
             write = TRUE)
load('five_params/forSTAT-5D-const-Q0-extended/ranges.Rdata')


### 10/22 Updates
combine_data(new_folder = 'to_Jake/MATTER+LBT_method-1/run_3_15-points/forSTAT-5D-15_extra_points/',
             old_folder = 'to_Jake/MATTER+LBT_method-1/combined_5D/',
             datasets = paste0(all_dsets, '.dat'),
             new_lhc = 'latin_hc_five_const_extra.txt',
             old_lhc = 'combined_lhc_5d.dat',
             new_ranges_file = 'ranges.Rdata',
             old_ranges_file = 'ranges.Rdata',
             write_folder = 'to_Jake/MATTER+LBT_method-1/combined_5D_1-3/',
             lhc_name = 'method_1_combined_1_3_lhc.txt',
             scaled_d_name = 'method_1_combined_1_3_scaled_d.txt',
             old_header = TRUE,
             new_header = FALSE,
             write = TRUE)


combine_data(new_folder = 'to_Jake/MATTER+LBT_method-2/run_2_12-points/forSTAT/',
            old_folder = 'to_Jake/MATTER+LBT_method-2/forSTAT-final-4D/',
            datasets = paste0(all_dsets, '.dat'),
            new_lhc = 'latin_hc_four_final_added.txt',
            old_lhc = 'latin_hc_four_final.txt',
            new_ranges_file = 'ranges.Rdata',
            old_ranges_file = 'ranges.Rdata',
            write_folder = 'to_Jake/MATTER+LBT_method-2/combined_4D_1-2/',
            lhc_name = 'method_2_combined_1_2_lhc.txt',
            scaled_d_name = 'method_2_combined_1_2_scaled_d.txt',
            old_header = FALSE,
            new_header = FALSE,
            write = TRUE)

##12/15 update
combine_data(new_folder = 'to_Jake/LBT/b_d_add_2/',
             old_folder = 'to_Jake/LBT/b_d_add_1/',
             datasets = paste0(all_dsets, '.dat'),
             new_lhc = 'b_d_add_2_lhc.dat',
             old_lhc = 'b_d_add_1_lhc.dat',
             new_ranges_file = 'b_d_add_2_ranges.Rdata',
             old_ranges_file = 'b_d_add_1_ranges.Rdata',
             write_folder = 'to_Jake/LBT/b_d_add_both/',
             lhc_name = 'fill_bd_lhc.dat',
             scaled_d_name = 'fill_bd_scaled_d.dat',
             old_header = FALSE,
             new_header = FALSE,
             write = TRUE)

##1/2 update
combine_data(new_folder = 'to_Jake/LBT/b_d_add_both/',
             old_folder = 'to_Jake/LBT/combine_all/',
             datasets = paste0(all_dsets, '.dat'),
             new_lhc = 'fill_bd_lhc.dat',
             old_lhc = 'combine_lhc_lbt_all.dat',
             new_ranges_file = 'b_d_add_2_ranges.Rdata',
             old_ranges_file = 'ranges.Rdata',
             write_folder = 'to_Jake/LBT/combine_all_bd_extra/',
             lhc_name = 'combine_all_bd_extra_lhc.dat',
             scaled_d_name = 'combine_all_bd_extra_scaled_d.dat',
             old_header = TRUE,
             new_header = FALSE,
             write = TRUE)

##2/3 update
combine_data(new_folder = 'LBT+MATTER_correction/LBT+MATTER-Method-2/forSTAT-4D-add-9pts/',
             old_folder = 'LBT+MATTER_correction/LBT+MATTER-Method-2/forSTAT-4D-60pts/',
             datasets = paste0(all_dsets, '.dat'),
             new_lhc = 'm2_add_2_1.txt',
             old_lhc = 'latin_hc_four_final.txt',
             new_ranges_file = 'ranges.Rdata',
             old_ranges_file = 'ranges.Rdata',
             write_folder = 'LBT+MATTER_correction/LBT+MATTER-Method-2/combine_m2_1-2/',
             lhc_name = 'combine_m2_lhc.dat',
             scaled_d_name = 'combine_m2_scaled_d.dat',
             old_header = FALSE,
             new_header = FALSE,
             write = TRUE)

View(read.table('LBT+MATTER_correction/LBT+MATTER-Method-2/combine_m2_1-2/AuAu200-cen-00-10.dat',header = FALSE))

View(read.table('LBT+MATTER_correction/LBT+MATTER-Method-2/forSTAT-4D-add-9pts/AuAu200-cen-00-10.dat',header = FALSE))

