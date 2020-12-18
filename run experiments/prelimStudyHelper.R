########
# prelimStudyHelper;
# run experiments multiple times
########
library(mclustcomp)
library(locfit)

makeMeasObj <- function(ARI, mclustObj, p, mse, z_d_sites) {
  return(list(ARI=ARI, mObj=mclustObj, purity=p, mse=mse, zds=z_d_sites))
}

ClusterPurity <- function(clusters, classes) {
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}

bundleStats <- function(pred_sites, og_sites, mse, z_d_sites, is_eBird=F){
  cmethods = c("jaccard", "mmm", "vi", "overlap", "nmi1")
  if(is_eBird){
    m <- mclustcomp(og_sites, og_sites, types=cmethods)
    clust_m_list <- makeMeasObj(0, m, 0, mse, 0)
  }else{
    ARI <- adjustedRandIndex(og_sites, pred_sites)
    
    m <- mclustcomp(pred_sites, og_sites, types=cmethods)
    p <- ClusterPurity(as.factor(pred_sites), as.factor(og_sites))
    
    clust_m_list <- makeMeasObj(ARI, m, p, mse, z_d_sites) 
  }
  
  return(clust_m_list)
}

measureDS <- function(){
  occ <- 0
  det <- 0
  ARI <- 0
  purity <- 0
  zds <- 0
  jacc <- 0
  mmm <- 0
  vi <- 0
  nmi <- 0
  overlap <- 0
  
  return(list(occ=occ, det=det, ARI=ARI, p=purity, jacc=jacc, mmm=mmm, vi=vi, nmi=nmi, overlap=overlap))
}

addToDS <- function(ds, stats_obj){
  ds$occ <- ds$occ + stats_obj$mse$occ
  ds$det <- ds$det + stats_obj$mse$det
  ds$ARI <- ds$ARI + stats_obj$ARI
  ds$p <- ds$p + stats_obj$purity
  ds$zds <- ds$zds + stats_obj$zds
  
  ds$jacc <- ds$jacc + stats_obj$mObj$scores[stats_obj$mObj$types == "jaccard"]
  ds$mmm <- ds$mmm + stats_obj$mObj$scores[stats_obj$mObj$types == "mmm"]
  ds$overlap <- ds$overlap + stats_obj$mObj$scores[stats_obj$mObj$types == "overlap"]
  ds$vi <- ds$vi + stats_obj$mObj$scores[stats_obj$mObj$types == "vi"]
  ds$nmi <- ds$nmi + stats_obj$mObj$scores[stats_obj$mObj$types == "nmi1"]
  
  return(ds)
}

divideDS <- function(ds, numIter){
  ds$occ <- ds$occ/numIter
  ds$det <- ds$det/numIter
  ds$ARI <- ds$ARI/numIter
  ds$p <- ds$p/numIter
  ds$zds <- ds$zds/numIter
  
  ds$jacc <- ds$jacc/numIter
  ds$mmm <- ds$mmm/numIter
  ds$overlap <- ds$overlap/numIter
  ds$vi <- ds$vi/numIter
  ds$nmi <- ds$nmi/numIter
  
  return(ds)
}

runIDExp <- function(kmsq_df, covObj, og_data, n_sites, geoClustAlpha=.8, occ_coef = TRUE_OCC_COEFF, det_coef = TRUE_DET_COEFF){
  
  occ_est <- 0
  i <- 0
  pvs_occ_obj <- NULL
  while(occ_est < 1000){
    i <- i + 1
    disp("iteration #", as.character(i))
    p_det <- 0
    while(p_det < .15){
      truth_df <- populateDF(kmsq_df, covObj$siteCovs, covObj$obsCovs, unique(kmsq_df$site), occ_coef, det_coef)
      og_data$species_observed_syn <- truth_df$species_observed_syn
      n_dets <- sum(truth_df$species_observed_syn)
      n_occ_sites <- sum(truth_df$occupied)
      
      p_det <- n_dets/nrow(truth_df)
      disp("number of detections: ", as.character(n_dets))
      disp("% detections: ", as.character(n_dets/nrow(truth_df)))
      # disp("% detected given its occupied: ", as.character(n_dets/n_occ_sites))
      disp("occupied probability: ", as.character(mean(truth_df$occupied_prob)))
    }
  
    clustGeoSites_obj_100000 <- runClustGeoID(truth_df, covObj, og_data, alpha=.8)
    # MAX_OBS = 10
    # clustGeoSites_obj_10 <- runClustGeoID(truth_df, covObj, og_data, alpha=.8)
    # MAX_OBS = 100000
    
    occ_est <- clustGeoSites_obj_100000$clust900$stats$mse$occ
      
    disp("# of zds for 900: ", as.character(clustGeoSites_obj_100000$clust900$stats$zds))
    disp("occ_est for same is: ", as.character(occ_est))
    
    if(occ_est < 10){
      pvs_occ_obj <- clustGeoSites_obj_100000
      pvs_n_occ_sites <- n_occ_sites
    }
  }
  # determine occupancy estimates
  disp("occupancy is huge; determine occupancy estimates")
  p_f <- clustGeoSites_obj_100000$clust900$pred_form
  
  clustGeoSites_obj_100000$clust900$df$occ_calc <- apply(clustGeoSites_obj_100000$clust900$df, 1, function(x) is_occ(x, p_f, covObj))
  
  return(list(c_100000=clustGeoSites_obj_100000, pvs=pvs_occ_obj, c_100000_occ=n_occ_sites, pvs_occ=pvs_n_occ_sites))
  # clustGeoSites_obj_100000$clust900$pred_form
}

is_occ <- function(row, pred_form, covObj){
  occ_sum <- pred_form@estimates["state"]@estimates[[1]]
  occ_cov_list <- covObj$siteCovs
  for(k in 1:length(occ_cov_list)){
    pred_form@estimates["state"]@estimates[[k+1]]
    occ_sum <- occ_sum + (pred_form@estimates["state"]@estimates[[k+1]] * as.numeric(row[[occ_cov_list[k]]]))
    occ_prob <- expit(occ_sum)
  }
  return(occ_prob)
}

runSingleExp <- function(truth_df, kmsq_df, covObj, og_data, n_sites, occ_coef = TRUE_OCC_COEFF, det_coef = TRUE_DET_COEFF){

  kmSqMSE <- calcOccMSE(truth_df, covObj, TRUE_OCC_COEFF, TRUE_DET_COEFF, syn_spec = TRUE, skip_closure = T)
  ret_df <- truth_df
  # feeding in the observations but not the sites
  og_data$species_observed_syn <- truth_df$species_observed_syn

  disp("clustGeo ... alpha & sites")
  clustGeoAlpha_obj <- runSingleClustGeoAlpha(truth_df, covObj, og_data, n_sites=599)
  clustGeoSites_obj <- runSingleClustGeoSites(truth_df, covObj, og_data, alpha=.8)
  
  disp("SKATER!!")
  SKATERSites_obj <- runSingleSKATERSites(truth_df, covObj, og_data)
  
  disp("dbscan ... eps & minPts")
  dbscan.eps <- runSingleDBSCAN.eps(truth_df, covObj, og_data, minPts = 6)
  dbscan.minPts <- runSingleDBSCAN.minPts(truth_df, covObj, og_data, eps = .001)
  
  disp("DBSC!")
  DBSC_df <- runDBSC(og_data, covObj)
  DBSC_MSE <- calcOccMSE(
    sites_df = DBSC_df, 
    covariate_object = covObj, 
    true_occ_coefficients = TRUE_OCC_COEFF, 
    true_det_coefficients = TRUE_DET_COEFF, 
    syn_spec = TRUE
  )
  
  li_DBSC_z_d_s <- calc_zero_det_sites(list(DBSC_MSE$checklists))
  
  DBSCStats <- bundleStats(DBSC_MSE$checklists$site, ret_df$site, DBSC_MSE$mse, li_DBSC_z_d_s[1], is_eBird = F)
  
  
  ##########
  # filter into sites
  ##########
  sites_ebird_filter <- filter_repeat_visits(
    og_data,
    min_obs = 2,
    max_obs = 10,
    annual_closure = TRUE,
    date_var = "formatted_date",
    site_vars = c("locality_id", "observer_id")
  )

  WETA_2017_4 <- roundLatLong(og_data, 4)
  sites_ebird_filter_ROUNDED <- filter_repeat_visits(
    WETA_2017_4,
    min_obs = MIN_OBS,
    max_obs = MAX_OBS,
    annual_closure = TRUE,
    date_var = "formatted_date",
    site_vars = c("rounded_locality_id")
  )

  ##########
  rounded_MSE <- calcOccMSE(
    sites_df = sites_ebird_filter_ROUNDED,
    covariate_object = covObj,
    true_occ_coefficients = TRUE_OCC_COEFF,
    true_det_coefficients = TRUE_DET_COEFF,
    syn_spec = TRUE
  )

  eBird_MSE <- calcOccMSE(
    sites_df = sites_ebird_filter,
    covariate_object = covObj,
    true_occ_coefficients = TRUE_OCC_COEFF,
    true_det_coefficients = TRUE_DET_COEFF,
    syn_spec = TRUE
  )
  
  li_eBird_z_d_s <- calc_zero_det_sites(list(eBird_MSE$checklists))
  li_round_z_d_s <- calc_zero_det_sites(list(rounded_MSE$checklists))
  li_kmSq_z_d_s <- calc_zero_det_sites(list(kmSqMSE$checklists))
  
  eBirdStats <- bundleStats(eBird_MSE$checklists$site, ret_df$site, eBird_MSE$mse, li_eBird_z_d_s, is_eBird=T)
  roundStats <- bundleStats(rounded_MSE$checklists$site, ret_df$site, rounded_MSE$mse, li_round_z_d_s)

  kmSqStats <- bundleStats(kmSqMSE$checklists, ret_df$site, kmSqMSE$mse, li_kmSq_z_d_s, is_eBird = T)
  
  
  return(list(kmSq=kmSqStats,
              eBird=eBirdStats, rounded=roundStats, 
              clust8=clustGeoAlpha_obj$clust8, clust5=clustGeoAlpha_obj$clust5, clust2=clustGeoAlpha_obj$clust2,
              clust900=clustGeoSites_obj$clust900, clust600=clustGeoSites_obj$clust600, clust300=clustGeoSites_obj$clust300,
              SKATER900=SKATERSites_obj$SKATER900, SKATER600=SKATERSites_obj$SKATER600, SKATER300=SKATERSites_obj$SKATER300,
              dbscan1=dbscan.eps$dbscan1, dbscan01=dbscan.eps$dbscan01, dbscan001=dbscan.eps$dbscan001,
              dbscan3=dbscan.minPts$dbscan3, dbscan6=dbscan.minPts$dbscan6, dbscan9=dbscan.minPts$dbscan9,
              DBSC=DBSCStats))
  
}


runMultExp <- function(kmsq_df, covObj, numIter, og_data, n_sites, geoClustAlpha=.8, occ_coef = TRUE_OCC_COEFF, det_coef = TRUE_DET_COEFF){
  eBird_DS <- measureDS()
  round_DS <- measureDS()
  
  clust8_DS <- measureDS()
  clust5_DS <- measureDS()
  clust2_DS <- measureDS()
  
  clust900_DS <- measureDS()
  clust600_DS <- measureDS()
  clust300_DS <- measureDS()
  
  SKATER900_DS <- measureDS()
  SKATER600_DS <- measureDS()
  SKATER300_DS <- measureDS()
  
  dbscan3_DS <- measureDS()
  dbscan6_DS <- measureDS()
  dbscan9_DS <- measureDS()
  
  dbscan1_DS <- measureDS()
  dbscan01_DS <- measureDS()
  dbscan001_DS <- measureDS()
  
  kmSq_DS <- measureDS()
  
  DBSC_DS <- measureDS()
  
  for(iter in 1:numIter){
    disp("beginning iteration #", iter)
    p_det <- 0
    while(p_det < .15){
      truth_df <- populateDF(kmsq_df, covObj$siteCovs, covObj$obsCovs, unique(kmsq_df$site), occ_coef, det_coef)
      
      n_dets <- sum(truth_df$species_observed_syn)
      n_occ_sites <- sum(truth_df$occupied)
      
      p_det <- n_dets/nrow(truth_df)
      disp("number of detections: ", as.character(n_dets))
      disp("% detections: ", as.character(n_dets/nrow(truth_df)))
      disp("% detected given its occupied: ", as.character(n_dets/n_occ_sites))
    }
    
    ret_obj <- runSingleExp(truth_df, kmsq_df, covObj = covObj, og_data, n_sites)
    MAX_OBS = 10
    ret_obj2 <- runSingleExp(truth_df, kmsq_df, covObj = covObj, og_data, n_sites)
    MAX_OBS = 100000
    
    ################
    # max_obs_100000
    ################
    eBird_DS <- addToDS(eBird_DS, ret_obj$eBird)
    round_DS <- addToDS(round_DS, ret_obj$rounded)
    
    clust8_DS <- addToDS(clust8_DS, ret_obj$clust8)
    clust5_DS <- addToDS(clust5_DS, ret_obj$clust5)
    clust2_DS <- addToDS(clust2_DS, ret_obj$clust2)
    
    clust900_DS <- addToDS(clust900_DS, ret_obj$clust900)
    clust600_DS <- addToDS(clust600_DS, ret_obj$clust600)
    clust300_DS <- addToDS(clust300_DS, ret_obj$clust300)
    
    SKATER900_DS <- addToDS(SKATER900_DS, ret_obj$SKATER900)
    SKATER600_DS <- addToDS(SKATER600_DS, ret_obj$SKATER600)
    SKATER300_DS <- addToDS(SKATER300_DS, ret_obj$SKATER300)
    
    dbscan3_DS <- addToDS(dbscan3_DS, ret_obj$dbscan3)
    dbscan6_DS <- addToDS(dbscan6_DS, ret_obj$dbscan6)
    dbscan9_DS <- addToDS(dbscan9_DS, ret_obj$dbscan9)
    
    dbscan1_DS <- addToDS(dbscan1_DS, ret_obj$dbscan1)
    dbscan01_DS <- addToDS(dbscan01_DS, ret_obj$dbscan01)
    dbscan001_DS <- addToDS(dbscan001_DS, ret_obj$dbscan001)
    DBSC_DS <- addToDS(DBSC_DS, ret_obj$DBSC)
    
    kmSq_DS <- addToDS(kmSq_DS, ret_obj$kmSq)
    ################
    
    ################
    # max_obs_10
    ################
    eBird_DS2 <- addToDS(eBird_DS2, ret_obj2$eBird)
    round_DS2 <- addToDS(round_DS2, ret_obj2$rounded)

    clust8_DS2 <- addToDS(clust8_DS2, ret_obj2$clust8)
    clust5_DS2 <- addToDS(clust5_DS2, ret_obj2$clust5)
    clust2_DS2 <- addToDS(clust2_DS2, ret_obj2$clust2)

    clust900_DS2 <- addToDS(clust900_DS2, ret_obj2$clust900)
    clust600_DS2 <- addToDS(clust600_DS2, ret_obj2$clust600)
    clust300_DS2 <- addToDS(clust300_DS2, ret_obj2$clust300)

    SKATER900_DS2 <- addToDS(SKATER900_DS2, ret_obj2$SKATER900)
    SKATER600_DS2 <- addToDS(SKATER600_DS2, ret_obj2$SKATER600)
    SKATER300_DS2 <- addToDS(SKATER300_DS2, ret_obj2$SKATER300)
    
    dbscan3_DS2 <- addToDS(dbscan3_DS2, ret_obj2$dbscan3)
    dbscan6_DS2 <- addToDS(dbscan6_DS2, ret_obj2$dbscan6)
    dbscan9_DS2 <- addToDS(dbscan9_DS2, ret_obj2$dbscan9)
    
    dbscan1_DS2 <- addToDS(dbscan1_DS2, ret_obj2$dbscan1)
    dbscan01_DS2 <- addToDS(dbscan01_DS2, ret_obj2$dbscan01)
    dbscan001_DS2 <- addToDS(dbscan001_DS2, ret_obj2$dbscan001)
    DBSC_DS2 <- addToDS(DBSC_DS2, ret_obj$DBSC)
    
    kmSq_DS2 <- addToDS(kmSq_DS2, ret_obj$kmSq)
  }  
  
  ################
  # max_obs_100000
  ################
  eBird_DS <- divideDS(eBird_DS, numIter)
  round_DS <- divideDS(round_DS, numIter)
  
  clust8_DS <- divideDS(clust8_DS, numIter)
  clust5_DS <- divideDS(clust5_DS, numIter)
  clust2_DS <- divideDS(clust2_DS, numIter)

  clust900_DS <- divideDS(clust900_DS, numIter)
  clust600_DS <- divideDS(clust600_DS, numIter)
  clust300_DS <- divideDS(clust300_DS, numIter)

  SKATER900_DS <- divideDS(SKATER900_DS, numIter)
  SKATER600_DS <- divideDS(SKATER600_DS, numIter)
  SKATER300_DS <- divideDS(SKATER300_DS, numIter)
  
  dbscan3_DS <- divideDS(dbscan3_DS, numIter)
  dbscan6_DS <- divideDS(dbscan6_DS, numIter)
  dbscan9_DS <- divideDS(dbscan9_DS, numIter)
  
  dbscan1_DS <- divideDS(dbscan1_DS, numIter)
  dbscan01_DS <- divideDS(dbscan01_DS, numIter)
  dbscan001_DS <- divideDS(dbscan001_DS, numIter)
  
  DBSC_DS <- divideDS(DBSC_DS, numIter)
  
  kmSq_DS <- divideDS(kmSq_DS, numIter)
  
  ################
  # max_obs_10
  ################
  eBird_DS <- divideDS(eBird_DS2, numIter)
  round_DS <- divideDS(round_DS2, numIter)
  
  clust8_DS2 <- divideDS(clust8_DS2, numIter)
  clust5_DS2 <- divideDS(clust5_DS2, numIter)
  clust2_DS2 <- divideDS(clust2_DS2, numIter)
  
  clust900_DS2 <- divideDS(clust900_DS2, numIter)
  clust600_DS2 <- divideDS(clust600_DS2, numIter)
  clust300_DS2 <- divideDS(clust300_DS2, numIter)
  
  SKATER900_DS2 <- divideDS(SKATER900_DS2, numIter)
  SKATER600_DS2 <- divideDS(SKATER600_DS2, numIter)
  SKATER300_DS2 <- divideDS(SKATER300_DS2, numIter)
  
  dbscan3_DS2 <- divideDS(dbscan3_DS2, numIter)
  dbscan6_DS2 <- divideDS(dbscan6_DS2, numIter)
  dbscan9_DS2 <- divideDS(dbscan9_DS2, numIter)
  
  dbscan1_DS2 <- divideDS(dbscan1_DS2, numIter)
  dbscan01_DS2 <- divideDS(dbscan01_DS2, numIter)
  dbscan001_DS2 <- divideDS(dbscan001_DS2, numIter)
  
  DBSC_DS2 <- divideDS(DBSC_DS2, numIter)
  
  kmSq_DS2 <- divideDS(kmSq_DS2, numIter)
  
  
  return(list(max_obs_100000=list(kmSq=kmSq_DS,
              eBird=eBird_DS, round=round_DS, 
              clust1=clust8_DS, clust5=clust5_DS, clust0=clust2_DS,
              clust900=clust900_DS, clust600=clust600_DS, clust300=clust300_DS,
              SKATER900=SKATER900_DS, SKATER600=SKATER600_DS, SKATER300=SKATER300_DS,
              dbscan3=dbscan3_DS, dbscan6=dbscan6_DS, dbscan9=dbscan9_DS,
              dbscan1=dbscan1_DS, dbscan01=dbscan01_DS, dbscan001=dbscan001_DS,
              DBSC=DBSC_DS),
              max_obs_10=list(kmSq=kmSq_DS2,
                   eBird=eBird_DS2, round=round_DS2, 
                   clust1=clust8_DS2, clust5=clust5_DS2, clust0=clust2_DS2,
                   clust900=clust900_DS2, clust600=clust600_DS2, clust300=clust300_DS2,
                   SKATER900=SKATER900_DS2, SKATER600=SKATER600_DS2, SKATER300=SKATER300_DS2,
                   dbscan3=dbscan3_DS2, dbscan6=dbscan6_DS2, dbscan9=dbscan9_DS2,
                   dbscan1=dbscan1_DS2, dbscan01=dbscan01_DS2, dbscan001=dbscan001_DS2,
                   DBSC=DBSC_DS2)))
}