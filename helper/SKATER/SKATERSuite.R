#####
# SKATER suite
#####


runSingleSKATERSites <- function(truth_df, covObj, og_data){       # can do num_sites_list <- c(900,600,300)
  ret_df <- truth_df
  truth_df <- subset(truth_df, select = -c(site))
  
  SKATER.MSE900 <- calcSKATER.MSE(
    og_data, 
    covObj, 
    numSites = 900, 
    enforce_false_p = FALSE
  )
  
  SKATER.MSE600 <- calcSKATER.MSE(
    og_data, 
    covObj, 
    numSites = 600, 
    enforce_false_p = FALSE
  )
  
  SKATER.MSE300 <- calcSKATER.MSE(
    og_data, 
    covObj, 
    numSites = 300, 
    enforce_false_p = FALSE
  )
  
  li_z_d_s <- calc_zero_det_sites(list(SKATER.MSE900$checklists, SKATER.MSE600$checklists, SKATER.MSE300$checklists))
  
  SKATER900Stats <- bundleStats(SKATER.MSE900$checklists$site, ret_df$site, SKATER.MSE900$mse, li_z_d_s[1])
  SKATER600Stats <- bundleStats(SKATER.MSE600$checklists$site, ret_df$site, SKATER.MSE600$mse, li_z_d_s[2])
  SKATER300Stats <- bundleStats(SKATER.MSE300$checklists$site, ret_df$site, SKATER.MSE300$mse, li_z_d_s[3])
  
  
  return(list(SKATER900=SKATER900Stats, SKATER600=SKATER600Stats, SKATER300=SKATER300Stats))
  
}
