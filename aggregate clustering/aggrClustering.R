########
# Aggregate Clustering based on this paper;
# https://dl.acm.org/doi/10.1145/1217299.1217303
########
# chunking based on MDD
########
library(sf)
library(sp)
library(raster)
setwd("/Users/MarkRoth/Documents/Oregon State/Research/eBird/occ and grouping checklists/occ-cluster/")
source("helper/DBSC/DBSCHelper.R")
source("helper/ClustGeo/clustGeoHelper.R")
source("helper/helpers.R")
source("helper/gen_sites.R")
source("aggregate clustering/aggClustHelper.R")
MDD <- 250    # meters
MIN_OBS <- 1
MAX_OBS <- 100000

# obj <- load.WETA()
# WETA_2017 <- obj[[1]]
# covObj <- obj[[2]]

TRUE_OCC_COEFF <- c(-.5, .85, 1, .2, -.5, -1)
TRUE_DET_COEFF <- c(1, 1, -.5, 1, -1, .5)

res_obj <- runExp()
# rm(res_obj)
# res_obj$round$mse
# res_obj$cg$mse
# res_obj$balls$mse
# res_obj$aggl$mse
# res_obj$base$mse


runExp <- function(){
  lst <- prepWETA.sites(WETA_2017, covObj, MDD, MIN_OBS, MAX_OBS)
  
  WETA_sites <- lst[[1]]
  og_data <- lst[[2]]
  proj_cent <- lst[[3]]
  truth_df <- lst[[4]]
  
  # generate sites!
  # tests <- list(rounded=list(4), eBird=NA, eBird_simple=T,
  # DBSC=T, clustGeo=list(c(.9, 850)))
  tests <- list(rounded=NA, eBird=NA, eBird_simple=T,DBSC=T,
                clustGeo=list(c(.9, 850), c(.5, 850), c(.1, 850)))
  WETA_sites <- appendSites(tests, WETA_sites, og_data)
  comb_df <- combineMethods(proj_cent, WETA_sites, og_data)
  comb_aggl_df <- combineMethodsAgg(proj_cent, WETA_sites, og_data)
  
  aggl_cons_df <- inner_join(comb_aggl_df, og_data, by = "vertex")
  
  # compare against single site dfs
  clustGeo_df_single <- clustGeoSites(alpha = .8, og_data, covObj, num_sites = 850)
  WETA_2017_i <- roundLatLong(og_data, 4)
  eBird_rounded_df <- filter_repeat_visits(
    WETA_2017_i,
    min_obs = 1,
    max_obs = 1000000,
    annual_closure = TRUE,
    date_var = "formatted_date",
    site_vars = c("rounded_locality_id")
  )
  
  
  # run this through an occupancy model
  consensus_mse <- calcOccMSE(comb_df, covObj, TRUE_OCC_COEFF, TRUE_DET_COEFF, syn_spec = T)
  aggl_consensus_mse <- calcOccMSE(aggl_cons_df, covObj, TRUE_OCC_COEFF, TRUE_DET_COEFF, syn_spec = T)
  
  cg1_mse <- calcOccMSE(clustGeo_df_single, covObj, TRUE_OCC_COEFF, TRUE_DET_COEFF, syn_spec = T)
  rounded4_mse <- calcOccMSE(eBird_rounded_df, covObj, TRUE_OCC_COEFF, TRUE_DET_COEFF, syn_spec = T)
  
  base_mse <- calcOccMSE(truth_df, covObj, TRUE_OCC_COEFF, TRUE_DET_COEFF, syn_spec = T)
  
  rounded4_mse$mse
  cg1_mse$mse
  consensus_mse$mse
  aggl_consensus_mse$mse
  base_mse$mse
  return(list(round=rounded4_mse, cg=cg1_mse, balls=consensus_mse, aggl=aggl_consensus_mse, base=base_mse))
}

