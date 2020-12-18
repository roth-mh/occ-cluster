###########
# populate syn data
###########
library(locfit)
library(mclust)

# occ_coefficients <- c(-.5, 0.4650489593471517, 0.25523159396930695, 0.8861590762863591, -0.20144247919650793, 0.3855035300948283)
# det_coefficients <- c(-1, -0.3605617900907053, 0.8502497538703003, -0.012596564563653434, 0.10894040227556012, 0.4170190133546442)

populateDF <- function(sites_df, occ_cov_list, det_cov_list, sites_list, occ_coefficients, det_coefficients){
  j<-0
  for(eBird_site in sites_list){
    j = j+1
    if(j %% 100 == 0){
      print(j)
    }
    checklists_at_site <- subset(sites_df[sites_df$site == eBird_site,])
    # choose random site to serve as the occ prob (and occ covariates) baseline
    
    # site_baseline <- checklists_at_site[which.max(checklists_at_site$time_observations_started),]
    
    occ_sum <- occ_coefficients[1]
    occ_cov_list <- covObj$siteCovs
    for(k in 1:length(occ_cov_list)){
      # find centroid of each occ cov and multiply it by the 
      avg_occ_cov <- mean(checklists_at_site[[occ_cov_list[k]]])
      occ_sum <- occ_sum + (occ_coefficients[k+1] * avg_occ_cov)
      checklists_at_site[[occ_cov_list[k]]] <- avg_occ_cov
    }
    
    det_cov_list <- covObj$obsCovs
    for(row in 1:nrow(checklists_at_site)){
      det_sum <- det_coefficients[1]
      
      for(det_cov_i in 1:length(det_cov_list)){
        det_cov <- det_cov_list[det_cov_i]
        det_sum <- det_sum + (det_coefficients[det_cov_i+1] * checklists_at_site[row,][[det_cov]])
      } 
      checklists_at_site[row,]$species_observed_syn <- rbinom(1, 1, expit(det_sum))
    }
    
    
    occ_prob <- expit(occ_sum)
    is_occupied <- rbinom(1, 1, occ_prob)
    checklists_at_site$occupied_prob <- occ_prob
    checklists_at_site$occupied <- is_occupied
    
    if(is_occupied == 0){
      checklists_at_site$species_observed_syn <- 0
    }
    
    
    if(j==1){
      clust_geo_df = checklists_at_site
    } else {
      clust_geo_df = rbind(clust_geo_df, checklists_at_site)
    }
    
  }
  # determine number of occupied sites
  occ_sites <- sqldf("SELECT occupied FROM clust_geo_df GROUP BY site")
  disp("number of occupied sites is ", as.character(sum(occ_sites)))
  return(clust_geo_df)
}
#########
# kmSq <- kmsq.MSE(WETA_df = WETA_2017, rad_m = 1000)
# # setup the ground truth
# truth_df <- populateDF(kmSq, covObj$siteCovs, covObj$obsCovs, unique(kmSq$site))
# ret_df <- truth_df
# ##########
# # setup to feed df into other clustering algs
# truth_df <- subset(truth_df, select = -c(site))
# ##########
# # filter into sites
# ##########
# sites_ebird_filter <- filter_repeat_visits(truth_df, min_obs = MIN_OBS, max_obs = MAX_OBS, annual_closure = TRUE, date_var = "formatted_date", site_vars = c("locality_id", "observer_id"))
# WETA_2017_4 <- roundLatLong(truth_df, 4)
# sites_ebird_filter_ROUNDED <- filter_repeat_visits(WETA_2017_4, min_obs = MIN_OBS, max_obs = MAX_OBS, annual_closure = TRUE, date_var = "formatted_date", site_vars = c("rounded_locality_id", "observer_id"))
# ##########
# 
# # nbs <- findNeighbors2(as.data.frame(truth_df), .001)
# # S.MSE <- calcSKATER.MSE2(truth_df, nbs, covObj, num_sites = 900)
# S.MSE <- calcSKATER.MSE(truth_df, covObj, numSites = 900, enforce_false_p = FALSE)
# clust.MSE <- calcClustGeoMSE(.8, truth_df, covObj, num_sites = 915, enforce_false_p = FALSE)
# 
# rounded_MSE <- calcOccMSE(sites_df = sites_ebird_filter_ROUNDED, covariate_object = covObj, true_occ_coefficients = TRUE_OCC_COEFF, true_det_coefficients = TRUE_DET_COEFF, syn_spec = TRUE)
# eBird_MSE <- calcOccMSE(sites_df = sites_ebird_filter, covariate_object = covObj, true_occ_coefficients = TRUE_OCC_COEFF, true_det_coefficients = TRUE_DET_COEFF, syn_spec = TRUE)
# 
# 
# eBird_MSE$mse
# rounded_MSE$mse
# clust.MSE$mse
# S.MSE$mse
# 
# 
# adjustedRandIndex(ret_df$site, clust.MSE$checklists$site)
# adjustedRandIndex(ret_df$site, S.MSE$checklists$site)
# # adjustedRandIndex(ret_df$site, kmSq.MSE$checklists$site)
# adjustedRandIndex(ret_df$site, eBird_MSE$checklists$site)
# adjustedRandIndex(ret_df$site, rounded_MSE$checklists$site)
# ##########
# 
# eBird_stats <- siteStatsDt(eBird_MSE$checklists, covObj = covObj)
# rounded_stats <- siteStatsDt(rounded_MSE$checklists, covObj = covObj)
# clust_stats <- siteStatsDt(clust.MSE$checklists, covObj = covObj)
# SKATER_stats <- siteStatsDt(S.MSE$checklists, covObj = covObj)
# 
# kmSq_stats <- siteStatsDt(kmSq.MSE$checklists, covObj = covObj)
# 
