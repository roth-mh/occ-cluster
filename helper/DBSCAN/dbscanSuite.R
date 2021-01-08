# dbscanSuite.R


runSingleDBSCAN.eps <- function(truth_df, covObj, og_data, minPts = 6){
  ret_df <- truth_df
  truth_df <- subset(truth_df, select = -c(site))
  
  dbscan.MSE1 <- dbscan.MSE(og_data, covObj, eps = .1, minPts = minPts)
  dbscan.MSE01 <- dbscan.MSE(og_data, covObj, eps = .01, minPts = minPts)
  dbscan.MSE001 <- dbscan.MSE(og_data, covObj, eps = .001, minPts = minPts)
  
  li_z_d_s <- calc_zero_det_sites(list(dbscan.MSE1$checklists, dbscan.MSE01$checklists, dbscan.MSE001$checklists))
  
  dbscan1Stats <- bundleStats(dbscan.MSE1$checklists$site, ret_df$site, dbscan.MSE1$mse, li_z_d_s[[1]])
  dbscan01Stats <- bundleStats(dbscan.MSE01$checklists$site, ret_df$site, dbscan.MSE01$mse, li_z_d_s[[2]])
  dbscan001Stats <- bundleStats(dbscan.MSE001$checklists$site, ret_df$site, dbscan.MSE001$mse, li_z_d_s[[3]])
  
  
  return(list(dbscan1=dbscan1Stats, dbscan01=dbscan01Stats, dbscan001=dbscan001Stats))
  
}

runSingleDBSCAN.minPts <- function(truth_df, covObj, og_data, eps = .001){
  ret_df <- truth_df
  truth_df <- subset(truth_df, select = -c(site))
  
  dbscan.MSE3 <- dbscan.MSE(og_data, covObj, eps = eps, minPts = 3)
  dbscan.MSE6 <- dbscan.MSE(og_data, covObj, eps = eps, minPts = 6)
  dbscan.MSE9 <- dbscan.MSE(og_data, covObj, eps = eps, minPts = 9)
  
  li_z_d_s <- calc_zero_det_sites(list(dbscan.MSE3$checklists, dbscan.MSE6$checklists, dbscan.MSE9$checklists))
  
  dbscan3Stats <- bundleStats(dbscan.MSE3$checklists$site, ret_df$site, dbscan.MSE3$mse, li_z_d_s[[1]])
  dbscan6Stats <- bundleStats(dbscan.MSE6$checklists$site, ret_df$site, dbscan.MSE6$mse, li_z_d_s[[2]])
  dbscan9Stats <- bundleStats(dbscan.MSE9$checklists$site, ret_df$site, dbscan.MSE9$mse, li_z_d_s[[3]])
  
  
  return(list(dbscan3=dbscan3Stats, dbscan6=dbscan6Stats, dbscan9=dbscan9Stats))
  
}
