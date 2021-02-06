########################
# ClustGeo helper
#   1. occ MSE calc
########################

#######
# Calculates the MSE of ClustGeo Grouping Algorithm
#######
calcClustGeoMSE <- function(alpha, checklists, covObj, occ_coeff, det_coeff, num_sites=NULL, ratio=NULL, filter = TRUE){
  if(is.null(num_sites) && is.null(ratio)){
    disp("NEED TO SPECIFY EITHER THE NUMBER OF SITES OR THE NUMBER OF CHECKLISTS/SITE")
    stopifnot(1 == 0)
  }
  
  # obsolete code
  if(filter){
    checklists <- filter_repeat_visits(
      checklists,
      min_obs = MIN_OBS, 
      max_obs = MAX_OBS,
      annual_closure = TRUE,
      date_var = "formatted_date",
      site_vars = c("locality_id")
    )
  }
  
  checklists_filtered <- subset(checklists, select = -c(site, n_observations, closure_id))  
  
  if(is.null(num_sites)){
    num_sites = round(nrow(checklists_filtered)/ratio)
  }
  
  WETA_env_data <- dist(subset(checklists_filtered, select = unlist(covObj$siteCovs)))
  WETA_geo_data <- dist(subset(checklists_filtered, select = c("latitude", "longitude")))
  
  tree <- hclustgeo(WETA_env_data, WETA_geo_data, alpha = alpha)
  part <- cutree(tree, num_sites)
  checklists_filtered$site <- part
  
  checklists <- data.frame(checklists_filtered)

  clustGeo_MSE <- calcOccMSE(sites_df = checklists,
                             covariate_object = covObj,
                             true_occ_coefficients = occ_coeff,
                             true_det_coefficients = det_coeff,
                             syn_spec = TRUE,
                             enforce_false_positives = TRUE
  )
  
  return(clustGeo_MSE)
}


#######
# helper function to link sites back 
# to all corresponding checklists
#######
find_site_clust <- function(row, checklists){
  checklists_filtered_df <- checklists
  i <- 1000000
  long <- row["longitude"]
  lat <- row["latitude"]
  site <- checklists_filtered_df[checklists_filtered_df$latitude == lat & checklists_filtered_df$longitdue == long,]$site
  if(nrow(site) == 0){
    site <- i
    i <- i + 1
  }
  return(as.character(site))
}
