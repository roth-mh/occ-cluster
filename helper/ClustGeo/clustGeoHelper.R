########################
# ClustGeo helper
#   1. occ MSE calc
########################

#######
# Calculates the MSE of ClustGeo Grouping Algorithm
#######
calcClustGeoMSE <- function(alpha, checklists, covObj, num_sites=NULL, ratio=NULL, filter = TRUE, uniq_loca = FALSE, enforce_false_p = TRUE, lam){
  if(is.null(num_sites) && is.null(ratio)){
    disp("NEED TO SPECIFY EITHER THE NUMBER OF SITES OR THE NUMBER OF CHECKLISTS/SITE")
    stopifnot(1 == 0)
  }
  
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
  
  if(uniq_loca){
    checklists_filtered <- subset(checklists, !duplicated(site))
    checklists_filtered <- subset(checklists_filtered, select = -c(site, n_observations, closure_id))
  } else {
    checklists_filtered <- subset(checklists, select = -c(site, n_observations, closure_id))  
  }
  
  
  if(is.null(num_sites)){
    num_sites = round(nrow(checklists_filtered)/ratio)
  }
  
  WETA_env_data <- dist(subset(checklists_filtered, select = unlist(covObj$occ_cov)))
  WETA_geo_data <- dist(subset(checklists_filtered, select = c("latitude", "longitude")))
  
  tree <- ?hclustgeo(WETA_env_data, WETA_geo_data, alpha = alpha)
  part <- cutree(tree, num_sites)
  checklists_filtered$site <- part
  
  checklists_filtered_df <- data.frame(checklists_filtered)
  
  if(uniq_loca){
    # add the sites back to the checklists with the same lat/long
    checklists$site <- apply(checklists, 1, find_site_clust, checklists_filtered_df)
  } else {
    checklists <- checklists_filtered_df
  }

  # clustGeo_MSE <- calcOccMSE(sites_df = checklists,
  #                            covariate_object = covObj,
  #                            true_occ_coefficients = TRUE_OCC_COEFF,
  #                            true_det_coefficients = TRUE_DET_COEFF,
  #                            syn_spec = TRUE,
  #                            enforce_false_positives = enforce_false_p
  # )
  
  clustGeo_MSE <- calcOccMSELambda(sites_df = checklists,
                             covariate_object = covObj,
                             true_occ_coefficients = TRUE_OCC_COEFF,
                             true_det_coefficients = TRUE_DET_COEFF,
                             syn_spec = TRUE,
                             enforce_false_positives = enforce_false_p,
                             lambda = lam
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
