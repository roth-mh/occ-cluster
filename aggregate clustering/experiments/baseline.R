####
# file to run baseline prelim experiments
####

library(sf)
library(sp)
library(raster)
library(dbscan)
library(cluster)
# setwd("/Users/MarkRoth/Documents/Oregon State/Research/eBird/occ and grouping checklists/occ-cluster/")
source("helper/DBSC/DBSCHelper.R")
source("helper/ClustGeo/clustGeoHelper.R")
source("helper/helpers.R")
source("helper/gen_sites.R")
source("helper/kmsq.R")
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

# only 1/alg for now
test_names <- list(
  "na1",
  "na2",
  "base",
  "clustGeo-.8-850",
  "eBird_simple"
  # "eBird",
  # "rounded",
  # "kmSq-1000",
  # "DBSC"
)

tests <- genTests(test_names)

results <- list(
  rounded=NA, 
  eBird=NA, 
  eBird_simple=NA,
  # kmSq=list(),
  kmSq=NA,
  DBSC=NA,
  GT=NA,
  noisy_gt=NA,
  # clustGeo=list(),
  clustGeo=NA,
  # agnes=list(),
  agnes=NA,
  kmeans=NA,
  base=NA
  # kmeans=list()
  # local=list()
)
numExps <- 20
for(i in 1:numExps){
  TRUE_OCC_COEFF <- runif(6, -1.5, 1.5)
  TRUE_DET_COEFF <- runif(6, -1.5, 1.5)
  truth_df <- populateDF(WETA_filtered, covObj$siteCovs, covObj$obsCovs, unique(WETA_filtered$site), TRUE_OCC_COEFF, TRUE_DET_COEFF)  
  
  WETA_2017$site <- -1
  WETA_2017$vertex <- seq(1:nrow(WETA_2017))
  
  og_data <- subset(WETA_2017, select = -c(site))
  
  og_data <- og_data[order(og_data$checklist_id),]
  truth_df <- truth_df[order(truth_df$checklist_id),]
  
  og_data$species_observed_syn <- truth_df$species_observed_syn
  stopifnot(og_data$checklist_id == truth_df$checklist_id)
  list_of_clusterings <- baselineExp(tests, og_data, covObj, truth_df)
  
  for(j in 1:length(list_of_clusterings)){
    clustering <- list_of_clusterings[[j]][[1]]
    clustering_name <- names(list_of_clusterings)[[j]]
    t_name <- strsplit(clustering_name, "-")[[1]]
    
    if(!is.na(clustering)){
      # if(names(clustering) %in% c("kmSq", "clustGeo")){
      #   # for each test in there
      #   disp("iterating over list")
      # } else {
      res <- clusterStats(clustering, truth_df, clustering_name, covObj, TRUE_OCC_COEFF, TRUE_DET_COEFF)
      # if(is.na(results[[t_name]])){
      if(i == 1){
        results[[t_name]] <- new("Clustering",
                 name=clustering_name,
                 occ.mse=res$mse$occ,
                 det.mse=res$mse$det,
                 checklist_df=clustering,
                 simil.GT=res$simil.GT)
      } else {
        results[[t_name]]@occ.mse <- results[[t_name]]@occ.mse + res$mse$occ
        results[[t_name]]@det.mse <- results[[t_name]]@det.mse + res$mse$det
        results[[t_name]]@simil.GT <- results[[t_name]]@simil.GT + res$simil.GT
      }
      # as.data.frame(res, row.names = test_name)
    }
  }
}

saveRDS(results, file = "clustering_results_obj.rds")

# res <- readRDS("aggregate clustering/results/clustering_results_obj.rds")
# for(object in res){
#   disp(object@name)
# }

