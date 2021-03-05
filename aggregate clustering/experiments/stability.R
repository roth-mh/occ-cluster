
library(sf)
library(sp)
library(raster)
library(dbscan)
library(cluster)
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

# hard to test this with clustGeo when we don't know howthe ballpark of unique locations/sites...
tests <- list(rounded=T, eBird=NA, 
                eBird_simple=T,
                DBSC=T,
                clustGeo=list(),
                GT=NA,
                noisy_gt=NA,
                agnes=list(),
                kmeans=list(),
                local=list()
                #ANOTHER HERE?,
)

test_names <- list("balls",
                   "agglom-fast", 
                   "rounded",
                   "eBird_simple",
                   "DBSC")

####################
#### PRE LOADING OBJ
####################
preLoad <- list()
preLoad <- readRDS("preLoad-3-3.rds")
TRUE_OCC_COEFF <- runif(6, -1.5, 1.5)
TRUE_DET_COEFF <- runif(6, -1.5, 1.5)

if(length(preLoad) == 0){
  lst <- prepProj.centers(WETA_2017, covObj, MDD, MIN_OBS, MAX_OBS)
  proj_cent <- lst[[1]]
  WETA_filtered <- lst[[2]]
} else {
  proj_cent <- preLoad[[1]]
  WETA_filtered <- preLoad[[2]]
}

truth_df <- populateDF(WETA_filtered, covObj$siteCovs, covObj$obsCovs, unique(WETA_filtered$site), TRUE_OCC_COEFF, TRUE_DET_COEFF)  
gen.new.df <- TRUE

numExps <- 10
sampleSizes <- c(500, 1000, 1500, 2000)
for(sampSize in sampleSizes){
  for(i in 1:numExps){
    
    WETA_2017$site <- -1
    WETA_2017$vertex <- seq(1:nrow(WETA_2017))
    
    og_data <- subset(WETA_2017, select = -c(site))
    og_data$species_observed_syn <- truth_df$species_observed_syn
    
    WETA_sites <- subset(og_data, select = c(checklist_id, vertex))
    
    
    samp <- sample(seq(1:nrow(truth_df)), sampSize)
    WETA_sites <- WETA_sites[samp,]

    og_data <- og_data[og_data$checklist_id %in% WETA_sites$checklist_id,]
    truth_df <- truth_df[truth_df$checklist_id %in% WETA_sites$checklist_id,]
    proj_cent <- proj_cent[proj_cent$checklist_id %in% WETA_sites$checklist_id,]
    
    res_obj <- runExp(tests, covObj, WETA_sites, og_data, truth_df, proj_cent)
    
    #####
    # similarity btwn aggl method
    # and GT partition
    #####
    comp_to_ground_truth_df_i <- makeSIMIL_TO_GT.DF(res_obj, test_names)
    
    #####
    # similarity btwn aggl method
    # and each input method
    #####
    simil_input_method_df_i <- makeCLUSTER_COMP.DF(res_obj, test_names)
    
    #####
    # mse values
    #####
    mse_df_i <- makeMSE.DF(res_obj, test_names)
    
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
    
    if(gen.new.df){
      TRUE_OCC_COEFF <- runif(6, -1.5, 1.5)
      TRUE_DET_COEFF <- runif(6, -1.5, 1.5)
      truth_df <- populateDF(WETA_filtered, covObj$siteCovs, covObj$obsCovs, unique(WETA_filtered$site), TRUE_OCC_COEFF, TRUE_DET_COEFF)      
    }
    
  }
  exp_type <- "stability" # MUST BE distinct OR overlap OR noisy OR pCorr non-spatial OR stability OR <empty>
  exp_name <- paste0("-", exp_type, "-ss-", sampSize)
  write.csv(comp_to_ground_truth_df, file=paste0("aggregate clustering/results/", exp_type, "/similarity_to_GT_Exp", exp_name, ".csv"))
  write.csv(simil_input_method_df, file=paste0("aggregate clustering/results/", exp_type, "/similarity_to_agglom_Exp", exp_name, ".csv"))
  write.csv(mse_df, file=paste0("aggregate clustering/results/", exp_type, "/mse", exp_name, ".csv"))
}




