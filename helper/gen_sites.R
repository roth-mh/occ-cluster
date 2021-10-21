###########
# populate syn data
###########
library(locfit)
library(mclust)

# populates the DF with occ/det values & probabilities
populateDF <- function(sites_df, occ_cov_list, det_cov_list, sites_list, occ_coefficients, det_coefficients){
  j<-0
  for(eBird_site in sites_list){
    j = j+1
    # if(j %% 100 == 0){
    #   print(j)
    # }
    checklists_at_site <- subset(sites_df[sites_df$site == eBird_site,])
    # choose random site to serve as the occ prob (and occ covariates) baseline
    
    # site_baseline <- checklists_at_site[which.max(checklists_at_site$time_observations_started),]
    
    occ_sum <- occ_coefficients[1]
    for(k in 1:length(occ_cov_list)){
      # find centroid of each occ cov and multiply it by the 
      avg_occ_cov <- mean(checklists_at_site[[occ_cov_list[k]]])
      occ_sum <- occ_sum + (occ_coefficients[k+1] * avg_occ_cov)
      checklists_at_site[[occ_cov_list[k]]] <- avg_occ_cov
    }
    
    for(row in 1:nrow(checklists_at_site)){
      det_sum <- det_coefficients[1]
      
      for(det_cov_i in 1:length(det_cov_list)){
        det_cov <- det_cov_list[det_cov_i]
        det_sum <- det_sum + (det_coefficients[det_cov_i+1] * checklists_at_site[row,][[det_cov]])
      } 
      det_prob <- expit(det_sum)
      checklists_at_site[row,]$det_prob <- det_prob
      checklists_at_site[row,]$species_observed_syn <- rbinom(1, 1, det_prob)
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


