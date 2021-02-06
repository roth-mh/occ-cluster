########
# Aggregate Clustering based on this paper;
# https://dl.acm.org/doi/10.1145/1217299.1217303
########
library(IntClust)


DBSC.sites <- runDBSC(WETA_2017, covObj)
DBSCSites <- DBSC.sites

length(unique(DBSC.sites$site))

# DBSC.mse <- calcOccMSE(DBSC.sites, covObj, TRUE_OCC_COEFF, TRUE_DET_COEFF, syn_spec = T)
# true.mse <- calcOccMSE(truth_df, covObj, TRUE_OCC_COEFF, TRUE_DET_COEFF, syn_spec = T)

# uniqLocations <- sqldf("SELECT site, latitude, longitude FROM DBSCSites GROUP BY latitude, longitude")
# visits_site <- sqldf("SELECT site, count(*) visits FROM uniqLocations GROUP BY site")

# at this point, do we need to compare all checklists at a site or only those that have different
# lat/lon points; can only 1 checklist represent all at a single lat/long?
# ...
# i think we can assume that all checklists at the same lat/long are in the site


# should do a thorough check on this... once syn species is figured out
# quick sanity check:
sites_ebird_filter <- filter_repeat_visits(
  WETA_2017,
  min_obs = 2,
  max_obs = 10,
  annual_closure = TRUE,
  date_var = "formatted_date",
  site_vars = c("locality_id", "observer_id")
)

eBird_MSE <- calcOccMSE(
  sites_df = sites_ebird_filter,
  covariate_object = covObj,
  true_occ_coefficients = TRUE_OCC_COEFF,
  true_det_coefficients = TRUE_DET_COEFF,
  syn_spec = TRUE
)

# a grouping only on same lat/long pairs
WETA_simple$species_observed_syn <- truth_df$species_observed_syn


simpleGrouped_MSE <- calcOccMSE(
  sites_df = WETA_simple,
  covariate_object = covObj,
  true_occ_coefficients = TRUE_OCC_COEFF,
  true_det_coefficients = TRUE_DET_COEFF,
  syn_spec = TRUE
)

eBird_MSE$mse
simpleGrouped_MSE$mse






########
# chunking based on MDD
########
library(sf)
library(sp)
library(raster)
setwd("/Users/MarkRoth/Documents/Oregon State/Research/eBird/occ and grouping checklists/occ-cluster/")
source("helper/DBSC/DBSCHelper.R")
source("helper/helpers.R")
source("helper/gen_sites.R")
source("aggregate clustering/aggClustHelper.R")
MDD <- 400    # meters
MIN_OBS <- 1
MAX_OBS <- 100000

TRUE_OCC_COEFF <- c(-.5, .85, 1, .2, -.5, -1)
TRUE_DET_COEFF <- c(1, 1, -.5, 1, -1, .5)

w_o <- load.WETA(all=T)
WETA_2017 <- w_o[[1]]
covObj <- w_o[[2]]

# WETA_filtered <- load.WETA_filtered(WETA_2017)
# truth_df <- populateDF(WETA_filtered, covObj$siteCovs, covObj$obsCovs, unique(WETA_filtered$site), TRUE_OCC_COEFF, TRUE_DET_COEFF)
# WETA_2017$species_observed_syn <- truth_df$species_observed_syn

WETA_uniq <- sqldf("SELECT * FROM WETA_2017 GROUP BY latitude, longitude")
WETA_uniq$site <- -1

ras.OR2 <- raster("~/Documents/Oregon State/Research/eBird/occ and grouping checklists/occ-cluster/site construction/input data/land cover data/OR_2020_fittedImage.tif")
prj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

################
# define the checklists as centroids of the buffer
# reproject centroids into appropriate projection for meters
# create buffer around each centroid of territory size
################
centroids <- st_as_sf(WETA_uniq, coords = c("longitude", "latitude"), crs=prj)
proj_cent <- st_transform(centroids, crs(ras.OR2))
buffers <- st_buffer(proj_cent, dist = MDD)

########
# for each buffer, find the checklists that overlap
# if all are unlabeled, create new cluster and add them
# if some are labeled (possibly some are unlabeled)
# union them
########
i <- nrow(WETA_uniq)
for(b in 1:nrow(buffers)){
  buff <- buffers[b,]
  over <- proj_cent[st_intersection(buff, proj_cent),]
  
  ########
  # plot all overlapping point for a given buffer
  ########
  # ggplot() +
  #   geom_sf(data = buff, pch = 4, col = "gray40") +
  #   geom_sf(data = over, col = "red") +
  #   theme_minimal()
  ########
  
  if(length(over) > 1){
    sites <- proj_cent[proj_cent$checklist_id %in% over$checklist_id,]$site
    # over$site
    # if sites are both -1: site <- i
    if(case1(sites)){
      proj_cent[proj_cent$checklist_id %in% over$checklist_id,]$site <- i
      i <- i + 1
    }
    # all of the same site
    else if(case2(sites)){
      
    }
    # some mix of sites ==> union them
    else {
      proj_cent <- chainSites(over, proj_cent, i)
      i <- i + 1
    }
    
  } else {
    # what is this condition? single point
    proj_cent[proj_cent$checklist_id == over$checklist_id,]$site <- i
    i <- i + 1
  }
}

# rename sites
j <- 1
for(s in unique(proj_cent$site)){
  proj_cent[proj_cent$site == s,]$site <- j
  j <- j + 1
}


hist(proj_cent$site, breaks=seq(0,length(unique(proj_cent$site))), main = sprintf("MDD: %s, # Visits/Site, all OR", MDD) )
dt <- as.data.frame(table(table(proj_cent$site)))
colnames(dt) <- c("num visits", "freq")
grid.table(dt)


# The goal is to group similar checklists
# So we pick out only certain checklists to compare
#   say we have a group of 10 checklists.
#     we make all pairwise comparisons and follow
#     the procedure detailed by the BALLS alg



