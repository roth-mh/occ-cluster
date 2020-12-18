######
# DBSCAN impl
######
library(dbscan)

dbscan.MSE <- function(WETA_df, covObj, eps, minPts){
  dbs <- WETA_df
  z <- subset(dbs, select = c(covObj$siteCovs, "latitude", "longitude"))
  dm<-data.matrix(z)
  
  dbscan <- dbscan::dbscan(dm, eps=eps, minPts = minPts)
  i <- max(dbscan$cluster) + 1
  for(j in 1:length(dbscan$cluster)){
    if(dbscan$cluster[j] == 0){
      dbscan$cluster[j] <- i
      i <- i + 1
    }
  }
  dbs$site <- dbscan$cluster
  
  z <- calcOccMSE(dbs, covObj, true_occ_coefficients = TRUE_OCC_COEFF, true_det_coefficients = TRUE_DET_COEFF, enforce_false_positives = F, syn_spec = T)
  return(z)
}

