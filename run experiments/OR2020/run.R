############
# file to run (sp. cl.) algs on eBird
# and compare results to silver standard
# (OR2020)
############
library(sf)
library(sp)
library(raster)
library(dbscan)
library(cluster)
library(tidyr)
setwd("/Users/MarkRoth/Documents/Oregon State/Research/eBird/occ and grouping checklists/occ-cluster/")
source("helper/DBSC/DBSCHelper.R")
source("helper/ClustGeo/clustGeoHelper.R")
source("helper/helpers.R")
source("helper/gen_sites.R")
source("helper/kmsq.R")
source("aggregate clustering/aggClustHelper.R")

OR.init <- stack("../TassledCapOR/OR_tasscap_summer_2011_cropped.tif")
names(OR.init) <- c("TCB", "TCG", "TCW", "TCA")


#######
# extract environmental features
# at checklist locations
#######
extractEnvFeat <- function(df, OR.tif){
  df.pts <- df
  coordinates(df.pts) <- c('longitude', 'latitude')
  crs(df.pts) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  
  pts_proj <- spTransform(df.pts, CRSobj=crs(OR.tif))
  
  final_df <- data.frame(
    base::scale(subset(df, select=
    -c(longitude, latitude, checklist_id, species_observed,sampling_event_identifier, 
      state_code, protocol_type, observation_date, latlong, observer_id, 
      scientific_name, locality_id, all_species_reported, Checklist_ID_Temp
    ))),
    subset(df, select=c(
      longitude, latitude, checklist_id, species_observed,sampling_event_identifier, 
      state_code, protocol_type, observation_date, latlong, observer_id, 
      scientific_name, locality_id, all_species_reported, Checklist_ID_Temp
    )),
    base::scale(raster::extract(OR.tif, pts_proj))
  ) 
  return(final_df)
}

#######
# calculates the occupancy model from a given dataset
# containing checklists and sites
#######
# 1. enforces closure
# 2. formats it w/r/t eBird data
# 3. runs through occupancy model
#######
calcOccModel <- function(df, occ_covs, det_covs, skip_closure=FALSE){
  sites_occ <- subset(df, !duplicated(site))$site
  # this (v v) function is synthetic species specific
  if(skip_closure){
    closed_df <- df
  } else {
    closed_df <- enforceClosure(df, site_covs, sites_occ)
  }
  
  umf_AUK <- auk::format_unmarked_occu(
    closed_df,
    site_id = "site",
    response = "species_observed",
    site_covs = occ_covs,
    obs_covs = det_covs
  )
  
  det_cov_str <- paste("", paste(det_covs, collapse="+"), sep=" ~ ")
  occ_cov_str <- paste("", paste(occ_covs, collapse="+"), sep=" ~ ")
  
  species_formula <- paste(det_cov_str, occ_cov_str, sep = " ")
  species_formula <- as.formula(species_formula)
  
  occ_um <- unmarked::formatWide(umf_AUK, type = "unmarkedFrameOccu")
  
  # og_syn_gen_form <- unmarked::occuPEN(formula = species_formula, occ_um, lambda = .1, pen.type = "Bayes")
  og_syn_gen_form <- unmarked::occu(formula = species_formula, occ_um)
  return(og_syn_gen_form)
}

########
# top level function to
# define and run experiment
########
runExp <- function(){
  CG1 <- "clustGeo-.8-850"
  # only 1/alg for now
  test_names <- list(
    "na",
    "na1",
    CG1,
    "eBird_simple",
    "eBird",
    "rounded-4",
    "kmSq-1000",
    "DBSC"
  )
  tests <- genTests(test_names)
  
  df <- ebd.no.or2020.2017_region
  
  baselineExp(df, tests, )
  
  prep.dataframe()
  
  df$occupied_prob <- -1
  df$det_prob <- -1
  df$formatted_date <- as.Date("01-01-2011")
  
  
  
  # (1) run checklists through algs, generate OM
  # (2) examine each OR2020 location.
  #       (i)  if occupied, compare against predicted occupied probability
  #       (ii) if not detected, compare against predicted occupied 
  #                 probability/detection probability (one should be low?)
}




