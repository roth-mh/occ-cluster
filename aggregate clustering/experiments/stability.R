
library(sf)
library(sp)
library(raster)
library(dbscan)
library(cluster)
# setwd("/Users/MarkRoth/Documents/Oregon State/Research/eBird/occ and grouping checklists/occ-cluster/")
source("helper/DBSC/DBSCHelper.R")
source("helper/ClustGeo/clustGeoHelper.R")
source("helper/helpers.R")
source("helper/gen_sites.R")
source("helper/kmsq.R")
source("aggregate clustering/aggClustHelper.R")

MDD <- 250    # meters
MIN_OBS <- 1
MAX_OBS <- 100000

obj <- load.WETA()
WETA_2017 <- obj[[1]]
covObj <- obj[[2]]


test_names <- list("balls",
                   "agglom-fast", 
                   "DBSC",
                   "eBird_simple",
                   "rounded-4")

tests <- genTests(test_names)

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
# TODO: CHANGE THIS SHIT
WETA_filtered_uniq <- sqldf("SELECT * from WETA_filtered GROUP BY latitude, longitude")

truth_df <- populateDF(WETA_filtered_uniq, covObj$siteCovs, covObj$obsCovs, unique(WETA_filtered_uniq$site), TRUE_OCC_COEFF, TRUE_DET_COEFF)  
gen.new.df <- TRUE

WETA_2017$site <- -1
WETA_2017$vertex <- seq(1:nrow(WETA_2017))

runStabilityExp(
  "boot", 
  truth_df = truth_df, 
  WETA_2017, 
  WETA_filtered_uniq, 
  proj_cent, 
  tests, 
  test_names, 
  deterministic = TRUE)






