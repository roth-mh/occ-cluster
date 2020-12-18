#clust Geo alpha investigation


runSingleClustGeoAlpha <- function(truth_df, covObj, og_data, n_sites){
  ret_df <- truth_df
  truth_df <- subset(truth_df, select = -c(site))
  
  clust.MSE8 <- calcClustGeoMSE(
    1, 
    og_data, 
    covObj, 
    num_sites = n_sites, 
    enforce_false_p = FALSE
  )
  
  disp("clust Geo 1 occ")
  disp(clust.MSE8$mse$occ)
  disp("clust Geo 1 det")
  disp(clust.MSE8$mse$det)
  
  clust.MSE5 <- calcClustGeoMSE(
    .5, 
    og_data, 
    covObj, 
    num_sites = n_sites, 
    enforce_false_p = FALSE
  )
  
  disp("clust Geo .5 occ")
  disp(clust.MSE5$mse$occ)
  disp("clust Geo .5 det")
  disp(clust.MSE5$mse$det)
  
  clust.MSE2 <- calcClustGeoMSE(
    .0, 
    og_data, 
    covObj, 
    num_sites = n_sites, 
    enforce_false_p = FALSE
  )
  
  li_z_d_s <- calc_zero_det_sites(list(clust.MSE8$checklists, clust.MSE5$checklists, clust.MSE2$checklists))
  
  clust8Stats <- bundleStats(clust.MSE8$checklists$site, ret_df$site, clust.MSE8$mse, li_z_d_s[1])
  clust5Stats <- bundleStats(clust.MSE5$checklists$site, ret_df$site, clust.MSE5$mse, li_z_d_s[2])
  clust2Stats <- bundleStats(clust.MSE2$checklists$site, ret_df$site, clust.MSE2$mse, li_z_d_s[3])
  
  
  return(list(clust8=clust8Stats, clust5=clust5Stats, clust2=clust2Stats))
  
}



runMultclustGeoAlphaExp <- function(kmsq_df, covObj, numIter, og_data, n_sites, geoClustAlpha=.8){
  clust8_DS <- measureDS()
  clust5_DS <- measureDS()
  clust2_DS <- measureDS()
  
  for(iter in 1:numIter){
    disp("beginning CLUSTGEO ALPHA iteration #", iter)
    
    ret_obj <- runSingleClustGeoAlpha(kmsq_df = kmsq_df, covObj = covObj, og_data, n_sites, geoClustAlpha)
    
    clust8_DS <- addToDS(clust8_DS, ret_obj$clust8)
    clust5_DS <- addToDS(clust5_DS, ret_obj$clust5)
    clust2_DS <- addToDS(clust2_DS, ret_obj$clust2)
  }  
  
  clust8_DS <- divideDS(clust8_DS, numIter)
  clust5_DS <- divideDS(clust5_DS, numIter)
  clust2_DS <- divideDS(clust2_DS, numIter)
  
  return(list(clust8=clust8_DS, clust5=clust5_DS, clust2=clust2_DS))
}




runSingleClustGeoSites <- function(truth_df, covObj, og_data, alpha){       # can do num_sites_list <- c(900,600,300)
  ret_df <- truth_df
  truth_df <- subset(truth_df, select = -c(site))
  
  clust.MSE900 <- calcClustGeoMSE(
    alpha, 
    og_data, 
    covObj, 
    num_sites = 900, 
    enforce_false_p = FALSE
  )
  
  clust.MSE600 <- calcClustGeoMSE(
    alpha, 
    og_data, 
    covObj, 
    num_sites = 600, 
    enforce_false_p = FALSE
  )
  
  clust.MSE300 <- calcClustGeoMSE(
    alpha, 
    og_data, 
    covObj, 
    num_sites = 300, 
    enforce_false_p = FALSE
  )
  
  li_z_d_s <- calc_zero_det_sites(list(clust.MSE900$checklists, clust.MSE600$checklists, clust.MSE300$checklists))
  
  clust900Stats <- bundleStats(clust.MSE900$checklists$site, ret_df$site, clust.MSE900$mse, li_z_d_s[1])
  clust600Stats <- bundleStats(clust.MSE600$checklists$site, ret_df$site, clust.MSE600$mse, li_z_d_s[2])
  clust300Stats <- bundleStats(clust.MSE300$checklists$site, ret_df$site, clust.MSE300$mse, li_z_d_s[3])
  
  
  return(list(clust900=clust900Stats, clust600=clust600Stats, clust300=clust300Stats))
  
}



runMultclustGeoAlphaExp <- function(kmsq_df, covObj, numIter, og_data, n_sites, geoClustAlpha=.8){
  clust8_DS <- measureDS()
  clust5_DS <- measureDS()
  clust2_DS <- measureDS()
  
  for(iter in 1:numIter){
    disp("beginning CLUSTGEO ALPHA iteration #", iter)
    
    ret_obj <- runSingleClustGeoAlpha(kmsq_df = kmsq_df, covObj = covObj, og_data, n_sites, geoClustAlpha)
    
    clust8_DS <- addToDS(clust8_DS, ret_obj$clust8)
    clust5_DS <- addToDS(clust5_DS, ret_obj$clust5)
    clust2_DS <- addToDS(clust2_DS, ret_obj$clust2)
  }  
  
  clust8_DS <- divideDS(clust8_DS, numIter)
  clust5_DS <- divideDS(clust5_DS, numIter)
  clust2_DS <- divideDS(clust2_DS, numIter)
  
  return(list(clust8=clust8_DS, clust5=clust5_DS, clust2=clust2_DS))
}






runClustGeoID <- function(truth_df, covObj, og_data, alpha){       # can do num_sites_list <- c(900,600,300)
  ret_df <- truth_df
  truth_df <- subset(truth_df, select = -c(site))
  
  clust.MSE900 <- calcClustGeoMSE(
    alpha, 
    og_data, 
    covObj, 
    num_sites = 900, 
    enforce_false_p = FALSE
  )
  
  clust.MSE600 <- calcClustGeoMSE(
    alpha, 
    og_data, 
    covObj, 
    num_sites = 600, 
    enforce_false_p = FALSE
  )
  
  clust.MSE300 <- calcClustGeoMSE(
    alpha, 
    og_data, 
    covObj, 
    num_sites = 300, 
    enforce_false_p = FALSE
  )
  
  li_z_d_s <- calc_zero_det_sites(list(clust.MSE900$checklists, clust.MSE600$checklists, clust.MSE300$checklists))
  
  clust900Stats <- bundleStats(clust.MSE900$checklists$site, ret_df$site, clust.MSE900$mse, li_z_d_s[1])
  clust600Stats <- bundleStats(clust.MSE600$checklists$site, ret_df$site, clust.MSE600$mse, li_z_d_s[2])
  clust300Stats <- bundleStats(clust.MSE300$checklists$site, ret_df$site, clust.MSE300$mse, li_z_d_s[3])
  
  
  return(list(clust900=list(stats=clust900Stats, df=clust.MSE900$checklists, pred_form=clust.MSE900$pred_form), clust600=list(stats=clust600Stats, df=clust.MSE600$checklists, pred_form=clust.MSE600$pred_form), clust300=list(stats=clust300Stats, df=clust.MSE300$checklists, pred_form=clust.MSE300$pred_form)))
  
}

