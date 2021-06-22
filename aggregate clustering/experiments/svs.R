# proof of concept that we can perform better than ground truth

library(sf)
library(sp)
library(raster)
library(dbscan)
library(cluster)
library(reshape2)
library(tidyr)
setwd("Documents/Oregon State/Research/eBird/occ and grouping checklists/occ-cluster/")
# source("helper/DBSC/DBSCHelper.R")
# source("helper/ClustGeo/clustGeoHelper.R")
source("helper/helpers.R")
source("helper/gen_sites.R")
# source("helper/kmsq.R")
source("aggregate clustering/aggClustHelper.R")

MDD <- 250    # meters
MIN_OBS <- 1
MAX_OBS <- 100000

obj <- load.WETA()
WETA_2017 <- obj[[1]]
covObj <- obj[[2]]

preLoad <- list()
# saveRDS(preLoad, file = "preLoad-3-3.rds")
preLoad <- readRDS("preLoad-3-3.rds")
# p_exp1 <- readRDS("pvsExp2-22.rds")

if(length(preLoad) == 0){
  lst <- prepProj.centers(WETA_2017, covObj, MDD, MIN_OBS, MAX_OBS)
  proj_cent <- lst[[1]]
  WETA_filtered <- lst[[2]]
} else {
  proj_cent <- preLoad[[1]]
  WETA_filtered <- preLoad[[2]]
}

# WETA filtered has checklists grouped into sites but has not yet enforced closure at the site

WETA_filtered.round <- WETA_filtered
WETA_filtered.round$latitude <- round(WETA_filtered.round$latitude, 5)
WETA_filtered.round$longitude <- round(WETA_filtered.round$longitude, 5)

WETA_filtered.round$latlong <- paste0(WETA_filtered.round$latitude, "_", WETA_filtered.round$longitude)

# disp(paste0("number of sites: ", n_sites))
n_svs <- 668

i <- 1
results <- list()
numExps <- 10
for(j in 1:numExps){
  disp("iteration #:", j)
  TRUE_OCC_COEFF <- runif(6, -1.5, 1.5)
  TRUE_DET_COEFF <- runif(6, -1.5, 1.5)
  truth_df <- populateDF(WETA_filtered, covObj$siteCovs, covObj$obsCovs, unique(WETA_filtered$site), TRUE_OCC_COEFF, TRUE_DET_COEFF)  
  
  n_svs <- 668
  for(exp_type in c("location_enforced", "random")){
    if(exp_type == "location_enforced"){
      svs.lim <- seq(670, 810, 10)
    } else {
      svs.lim <- append(seq(670, 810, 10), seq(846, nrow(truth_df), 100))
    }
    # defining the number of svs 
    svs.lim <- c(670,680,690)
    
    # iterating over each item in the list of svs
    for(n.svs.goal in svs.lim){
      disp("svs.goal is: ", n.svs.goal, " with experiment: ", exp_type)
      
      # resetting w_f to be the rounded filter and sorting w_f & truth
      w_f <- WETA_filtered.round
      w_f <- w_f[order(w_f$checklist_id),]
      truth_df <- truth_df[order(truth_df$checklist_id),]
      stopifnot(w_f$checklist_id == truth_df$checklist_id)
      
      # feeding in the species observations from the simulated dataset, truth_df
      w_f$species_observed_syn <- truth_df$species_observed_syn
      
      if(exp_type == "location_enforced" && n.svs.goal == 810){
        # results[[exp]][[as.character(n.svs.goal)]] <- results[[exp]][[as.character(930)]]
        n.svs.goal <- 809

      }
      
      WETA_sites <- sqldf("SELECT site, count(site) visits from w_f GROUP BY site")
      WETA_mvs <- WETA_sites[WETA_sites$visits > 1,]
      WETA_svs <- WETA_sites[WETA_sites$visits == 1,]
      
      n_svs <- nrow(WETA_svs)

      # while the goal is not met
      while(n_svs < n.svs.goal){
        
        if(exp_type == "location_enforced"){
          w_filtered <- filterOutSVS(w_f)
        }
        
        # ensuring we are changing a site with multiple visits to a single visit site & mvs (or svs if mvs had 2 checklists)
        n_sites <- length(unique(w_f$site))
        
        # Q1: does it matter if checklists are at different locations? 
        #         (w_f should have more svs because more mvs have only 2 visits? -- not necessarily true, unsure about the stats)
        # Q2: more svs or more mvs? (enforce svs?)
        if(exp_type == "location_enforced"){
          
          mvs <- w_filtered[w_filtered$site %in% WETA_mvs$site,]
          
          latlong_df <- sqldf("SELECT latlong, count(*) count FROM mvs GROUP BY latlong")
          uniq_latlong <- latlong_df[latlong_df$count == 1,]
          if(nrow(uniq_latlong) == 0){
            disp("we have finished")
          }
          to_change <- sample(uniq_latlong$latlong, 1)
          
          w_f[w_f$latlong == to_change,]$site <- as.character(paste0("s", i))
          
        } else {
          #######
          # not enforcing different location 
          #######
          mvs <- sample(unique(WETA_mvs$site), 1)
          mvs_checklists <- w_f[w_f$site == mvs,]
          row <- sample(nrow(mvs_checklists), 1)
          chklst <- mvs_checklists[row,]
          w_f[w_f$checklist_id == chklst$checklist_id,]$site <- as.character(paste0("s", i))
          
        }
        
        WETA_sites <- sqldf("SELECT site, count(site) visits from w_f GROUP BY site")
        WETA_svs <- WETA_sites[WETA_sites$visits == 1,]
        WETA_mvs <- WETA_sites[WETA_sites$visits > 1,]
        n_svs <- nrow(WETA_svs)
        i <- i + 1
        
      }
      
      alt.mse <- calcOccMSE(
        w_f,
        covObj,
        true_occ_coefficients = TRUE_OCC_COEFF,
        true_det_coefficients = TRUE_DET_COEFF,
        syn_spec = T)
      
      if(j == 1){
        results[[exp_type]][[as.character(n.svs.goal)]]$occ <- list(alt.mse$mse$occ)
        results[[exp_type]][[as.character(n.svs.goal)]]$det <- list(alt.mse$mse$det)
      } else {
        results[[exp_type]][[as.character(n.svs.goal)]]$occ <- append(results[[exp_type]][[as.character(n.svs.goal)]]$occ, alt.mse$mse$occ)
        results[[exp_type]][[as.character(n.svs.goal)]]$det <- append(results[[exp_type]][[as.character(n.svs.goal)]]$det, alt.mse$mse$det)
      }
    }
  }
  
  truth.mse <- calcOccMSE(
    truth_df,
    covObj,
    true_occ_coefficients = TRUE_OCC_COEFF,
    true_det_coefficients = TRUE_DET_COEFF,
    syn_spec = T)
  if(j == 1){
    results[["truth"]][[as.character(n.svs.goal)]]$occ <- list(truth.mse$mse$occ)
    results[["truth"]][[as.character(n.svs.goal)]]$det <- list(truth.mse$mse$det)
  } else {
    results[["truth"]][[as.character(n.svs.goal)]]$occ <- append(results[["truth"]][[as.character(n.svs.goal)]]$occ, truth.mse$mse$occ)
    results[["truth"]][[as.character(n.svs.goal)]]$det <- append(results[["truth"]][[as.character(n.svs.goal)]]$det, truth.mse$mse$det)
  }
  
}


# saveRDS(results, "SVSresults.RDS")
# results <- readRDS("aggregate clustering/results/cluster-results/SVSresults.RDS")

# results$truth$`2146`
cpy_results <- results
# results <- cpy_results
for(i in seq(670, 2146, 10)){
  results[["truth"]][[as.character(i)]] <- results[["truth"]][[1]]
}

res_mean <- list()
for(exp_i in names(results)){
  for(n_sites in names(results[[exp_i]])){
    res_mean[[as.character(exp_i)]][[as.character(n_sites)]]$occ <- paste0(mean(as.numeric(results[[exp_i]][[n_sites]]$occ)), "#", std(as.numeric(results[[exp_i]][[n_sites]]$occ)))
    # res_mean[[as.character(exp_i)]][[as.character(n_sites)]]$occ_std <- std(as.numeric(results[[exp_i]][[n_sites]]$occ))
#     
    res_mean[[as.character(exp_i)]][[as.character(n_sites)]]$det <- paste0(mean(as.numeric(results[[exp_i]][[n_sites]]$det)), "#", std(as.numeric(results[[exp_i]][[n_sites]]$det)))
    # res_mean[[as.character(exp_i)]][[as.character(n_sites)]]$det_std <- std(as.numeric(results[[exp_i]][[n_sites]]$det))
  }
}
# 
# 
melt.results <- reshape2::melt(res_mean)
# 
occ_res <- melt.results[melt.results$L3 == "occ",]
det_res <- melt.results[melt.results$L3 == "det",]
# 

occ_res <- occ_res %>%
  separate(value, c("mean", "std"), "#")
# 
det_res <- det_res %>%
  separate(value, c("mean", "std"), "#")
# 
# 
colnames(occ_res) <- c("mean", "std", "d/o", "num svs", "type")
colnames(det_res) <- c("mean", "std", "d/o", "num svs", "type")
# 
occ_res$mean <- as.numeric(occ_res$mean)
occ_res$std <- as.numeric(occ_res$std)
occ_res$`num svs` <- as.numeric(occ_res$`num svs`)
# 
det_res$mean <- as.numeric(det_res$mean)
det_res$std <- as.numeric(det_res$std)
det_res$`num svs` <- as.numeric(det_res$`num svs`)
# 
library(ggplot2)
ggplot(occ_res,
       aes(x = `num svs`,
           y = `mean`,
           color = type, group = type)) +
  geom_line() +
  geom_errorbar(aes(ymin=mean - std, ymax=mean+std )) +
  ggtitle("num single visit sites vs occupancy estimate 10 runs")

occ_res.zoom <- occ_res[occ_res$`num svs` < 810,]
ggplot(occ_res.zoom,
       aes(x = `num svs`,
       y = `mean`,
       color = type, group = type)) +
  geom_line() +
  geom_errorbar(aes(ymin=mean - std, ymax=mean+std )) +
  ggtitle("num single visit sites vs occupancy estimate 10 runs")

ggplot(det_res,
       aes(x = `num svs`,
           y = `mean`,
           color = type, group = type)) +
  geom_line() +
  geom_errorbar(aes(ymin=mean - std, ymax=mean+std )) +
  ggtitle("num single visit sites vs detection estimate 10 runs")

det_res.zoom <- det_res[det_res$`num svs` < 810,]
ggplot(det_res.zoom,
       aes(x = `num svs`,
           y = `mean`,
           color = type, group = type)) +
  geom_line() +
  geom_errorbar(aes(ymin=mean - std, ymax=mean+std )) +
  ggtitle("num single visit sites vs occupancy estimate 10 runs")






