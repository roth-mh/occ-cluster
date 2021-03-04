####
# file to test contribution of each alg
####


library(sf)
library(sp)
library(raster)
# library(dbscan)
setwd("/Users/MarkRoth/Documents/Oregon State/Research/eBird/occ and grouping checklists/occ-cluster/")
source("helper/DBSC/DBSCHelper.R")
source("helper/ClustGeo/clustGeoHelper.R")
source("helper/helpers.R")
source("helper/gen_sites.R")
source("aggregate clustering/aggClustHelper.R")

MDD <- 250    # meters
MIN_OBS <- 1
MAX_OBS <- 100000

obj <- load.WETA()
WETA_2017 <- obj[[1]]
covObj <- obj[[2]]

TRUE_OCC_COEFF <- c(-.5, .85, 1, .2, -.5, -1)
TRUE_DET_COEFF <- c(1, 1, -.5, 1, -1, .5)

p_exp <- readRDS(file="pvsExp2-22.rds")

# num runs / val
numRuns <- 10
for(p in c(.1,.3,.5,.7,.9)){
  # TO EXCLUDE, INPUT MUST BE NA or empty list for CG
  tests <- list(rounded=NA, 
                eBird=NA, 
                eBird_simple=T,
                DBSC=NA,
                clustGeo=list(),
                noisy_gt=p)
  
  test_names <- list("balls", 
                     "agglomerative", 
                     "ebdSimple",
                     paste0("GT-", as.character(p)))
  
  for(i in 1:numRuns){
    disp(paste0("running with GT noise == ", as.character(p), 
                " and iteration #", as.character(i)))
    
    res_obj2 <- runExp(tests, covObj, pvsExp = p_exp)
    
    # saving a pvsExp with rel-low mse values
    # saveRDS(res_obj2$pvsExp, file="pvsExp2-22.rds")
    
    #####
    # similarity btwn aggl method
    # and GT partition
    #####
    comp_df_i <- makeSIMIL_TO_GT.DF(res_obj2)
    
    #####
    # similarity btwn aggl method
    # and each input method
    #####
    simil_df_i <- makeCLUSTER_COMP.DF(res_obj2)
    
    #####
    # mse values
    #####
    mse_df_i <- makeMSE.DF(res_obj2)
    
    if(i != 1){
      temp <- cbind(comp_to_ground_truth_df, comp_df_i)
      comp_to_ground_truth_df <- sapply(unique(colnames(temp)), function(x) rowSums(temp[, colnames(temp) == x, drop = FALSE]))  
      
      temp <- cbind(simil_input_method_df, simil_df_i)
      simil_input_method_df <- sapply(unique(colnames(temp)), function(x) rowSums(temp[, colnames(temp) == x, drop = FALSE]))  
      
      temp <- cbind(mse_df, mse_df_i)
      mse_df <- sapply(unique(colnames(temp)), function(x) rowSums(temp[, colnames(temp) == x, drop = FALSE]))  
    } else {
      comp_to_ground_truth_df <- comp_df_i
      simil_input_method_df <- simil_df_i
      mse_df <- mse_df_i
    }
    if(p == 0){
      break
    }
  }
  
  if(p != 0){
    comp_to_ground_truth_df <- comp_to_ground_truth_df/10
    simil_input_method_df <- simil_input_method_df/10
    mse_df <- mse_df/10
  }
  
  exp_type <- "noisy" # MUST BE distinct OR overlap OR noisy OR none
  exp_num <- p
  exp_name <- paste0("-", exp_type, "-", exp_num)
  write.csv(comp_to_ground_truth_df, file=paste0("aggregate clustering/results/", exp_type, "/similarity_to_GT_Exp", exp_name, ".csv"))
  write.csv(simil_input_method_df, file=paste0("aggregate clustering/results/", exp_type, "/similarity_to_agglom_Exp", exp_name, ".csv"))
  write.csv(mse_df, file=paste0("aggregate clustering/results/", exp_type, "/mse", exp_name, ".csv"))
}

# for(tn in 1:length(test_names)){
#   disp(as.character(res_obj2$mse.list[[tn]]$mse))
# }

