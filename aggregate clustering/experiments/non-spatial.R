####
# file to test incorp. of non-spatial clustering algs
# (also test the overlap vs. new information balance)
####


library(sf)
library(sp)
library(raster)
library(dbscan)
library(cluster)
setwd("/Users/MarkRoth/Documents/Oregon State/Research/eBird/occ and grouping checklists/occ-cluster/")
source("helper/DBSC/DBSCHelper.R")
source("helper/ClustGeo/clustGeoHelper.R")
source("helper/helpers.R")
source("helper/gen_sites.R")
source("aggregate clustering/aggClustHelper.R")

MDD <- 250    # meters
MIN_OBS <- 1
MAX_OBS <- 100000

obj <- load.WETA()
WETA_2017 <- obj[[1]]
covObj <- obj[[2]]


test_names.A <- list("balls",
                     "agglom-updated",
                     "clustGeo-.8-850",
                     "DBSC",
                     "eBird_simple",
                     "all_svs")

test_names.B <- list("balls", 
                     "agglom-updated", 
                     "agnes-850",
                     "clustGeo-.8-850",
                     "DBSC",
                     "eBird_simple")

# xyz_names <- list("balls", 
#                      "agglom-updated", 
#                      "agnes-150",
#                      "clustGeo-.8-150",
#                      "DBSC",
#                      "eBird_simple",
#                      "rounded-4")
# tests.xyz <- genTests(xyz_names)


tests.A <- genTests(test_names.A)
tests.B <- genTests(test_names.B)


####################
#### PRE LOADING OBJ
####################
preLoad <- list()
# saveRDS(preLoad, file = "preLoad-3-3.rds")
preLoad <- readRDS("preLoad-3-3.rds")
# p_exp1 <- readRDS("pvsExp2-22.rds")
TRUE_OCC_COEFF <- runif(6, -1.5, 1.5)
TRUE_DET_COEFF <- runif(6, -1.5, 1.5)


if(length(preLoad) == 0){
  lst <- prepProj.centers(WETA_2017, covObj, MDD, MIN_OBS, MAX_OBS)
  proj_cent <- lst[[1]]
  WETA_filtered <- lst[[2]]
} else {
  proj_cent <- preLoad[[1]]
  WETA_filtered <- preLoad[[2]]
}

WETA_filtered_uniq <- sqldf("SELECT * from WETA_filtered GROUP BY latitude, longitude")
truth_df <- populateDF(WETA_filtered_uniq, covObj$siteCovs, covObj$obsCovs, unique(WETA_filtered_uniq$site), TRUE_OCC_COEFF, TRUE_DET_COEFF)  
gen.new.df <- TRUE
local_test <- FALSE

results.A <- list()
results.B <- list()

numExps <- 20
aggl_df.A <- NA
balls_df.A <- NA
aggl_df.B <- NA
balls_df.B <- NA
for(i in 1:numExps){
  
  WETA_filtered_uniq$site <- -1
  WETA_filtered_uniq$vertex <- seq(1:nrow(WETA_filtered_uniq))
  
  og_data <- subset(WETA_filtered_uniq, select = -c(site))
  
  og_data <- og_data[order(og_data$checklist_id),]
  truth_df <- truth_df[order(truth_df$checklist_id),]
  
  og_data$species_observed_syn <- truth_df$species_observed_syn
  
  WETA_sites <- subset(og_data, select = c(checklist_id, vertex))
  
  if(local_test){
    samp <- sample(seq(1:nrow(truth_df)), 600)
    W_sites <- WETA_sites[samp,]

    og_d <- og_data[og_data$checklist_id %in% W_sites$checklist_id,]
    t_df <- truth_df[truth_df$checklist_id %in% W_sites$checklist_id,]
    p_cent <- proj_cent[proj_cent$checklist_id %in% W_sites$checklist_id,]
  } else {
    W_sites <- WETA_sites
    og_d <- og_data
    t_df <- truth_df
    p_cent <- proj_cent
  }
  
  stopifnot(og_d$checklist_id == t_df$checklist_id)             
  res_obj.A <- runExp(tests.A, covObj, W_sites, og_d, t_df, p_cent, comb_df = balls_df.A, comb_aggl_df_fast = aggl_df.A)
  res_obj.B <- runExp(tests.B, covObj, W_sites, og_d, t_df, p_cent, comb_df = balls_df.B, comb_aggl_df_fast = aggl_df.B)
  
  aggl_df.A <- res_obj.A$aggl_df
  balls_df.A <- res_obj.A$balls_df
  aggl_df.B <- res_obj.B$aggl_df
  balls_df.B <- res_obj.B$balls_df

  for(j in 1:length(res_obj.A$sites_list)){
    clustering <- res_obj.A$sites_list[j][[1]]
    t_name <- names(res_obj.A$sites_list)[[j]]
    
    if(!is.na(clustering)){
      res <- clusterStats(clustering, t_df, t_name, covObj, TRUE_OCC_COEFF, TRUE_DET_COEFF)
      if(i == 1){
        results.A[[t_name]] <- new("Clustering",
                                 name=t_name,
                                 occ.mse=res$mse$occ,
                                 det.mse=res$mse$det,
                                 svs=res$svs,
                                 acc.svs=res$acc.svs,
                                 checklist_df=clustering,
                                 simil.GT=res$simil.GT)
      } else {
        results.A[[t_name]]@occ.mse <- results.A[[t_name]]@occ.mse + res$mse$occ
        results.A[[t_name]]@det.mse <- results.A[[t_name]]@det.mse + res$mse$det
        results.A[[t_name]]@svs <- results.A[[t_name]]@svs + res$svs
        results.A[[t_name]]@acc.svs <- results.A[[t_name]]@acc.svs + res$acc.svs
        results.A[[t_name]]@simil.GT <- results.A[[t_name]]@simil.GT + res$simil.GT
      }
    }
  }
  
  for(j in 1:length(res_obj.B$sites_list)){
    clustering <- res_obj.B$sites_list[j][[1]]
    t_name <- names(res_obj.B$sites_list)[[j]]
    
    if(!is.na(clustering)){
      res <- clusterStats(clustering, t_df, t_name, covObj, TRUE_OCC_COEFF, TRUE_DET_COEFF)
      # if(is.na(results[[t_name]])){
      if(i == 1){
        results.B[[t_name]] <- new("Clustering",
                                 name=t_name,
                                 occ.mse=res$mse$occ,
                                 det.mse=res$mse$det,
                                 svs=res$svs,
                                 acc.svs=res$acc.svs,
                                 checklist_df=clustering,
                                 simil.GT=res$simil.GT)
      } else {
        results.B[[t_name]]@occ.mse <- results.B[[t_name]]@occ.mse + res$mse$occ
        results.B[[t_name]]@det.mse <- results.B[[t_name]]@det.mse + res$mse$det
        results.B[[t_name]]@svs <- results.B[[t_name]]@svs + res$svs
        results.B[[t_name]]@acc.svs <- results.B[[t_name]]@acc.svs + res$acc.svs
        results.B[[t_name]]@simil.GT <- results.B[[t_name]]@simil.GT + res$simil.GT
      }
    }
  }

  if(gen.new.df){
    TRUE_OCC_COEFF <- runif(6, -1.5, 1.5)
    TRUE_DET_COEFF <- runif(6, -1.5, 1.5)
    truth_df <- populateDF(WETA_filtered_uniq, covObj$siteCovs, covObj$obsCovs, unique(WETA_filtered_uniq$site), TRUE_OCC_COEFF, TRUE_DET_COEFF)      
  }
  
}

z.A <- combineIntoDF(results.A, numExps)
write.csv(z.A, "non-spatial-20exps-A.csv")

z.B <- combineIntoDF(results.B, numExps)
write.csv(z.B, "non-spatial-20exps-B.csv")