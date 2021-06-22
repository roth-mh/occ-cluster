####
# file to run baseline prelim experiments
####

library(sf)
library(sp)
library(raster)
library(dbscan)
library(cluster)
library(tidyr)
setwd("/Users/MarkRoth/Documents/Oregon State/Research/eBird/occ and grouping checklists/occ-cluster/")
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
  "eBird_simple",
  "eBird",
  # "eBird_upper",
  # "eBird_lower"
  "rounded-4",
  "kmSq-1000",
  "DBSC"
)

# rounded affects about a ~50m^2 region

tests <- genTests(test_names)

# results <- list(
#   rounded=NA, 
#   eBird=NA, 
#   eBird_simple=NA,
#   # kmSq=list(),
#   kmSq=NA,
#   DBSC=NA,
#   GT=NA,
#   noisy_gt=NA,
#   # clustGeo=list(),
#   clustGeo=NA,
#   # agnes=list(),
#   agnes=NA,
#   kmeans=NA,
#   base=NA
#   # kmeans=list()
#   # local=list()
# )
# WETA_filtered_uniq <- sqldf("SELECT * from WETA_filtered GROUP BY latitude, longitude")

TRUE_OCC_COEFF <- runif(6, -1.5, 1.5)
TRUE_DET_COEFF <- runif(6, -1.5, 1.5)

results <- list()
numExps <- 2
for(exp_i in 1:numExps){
  disp("iteration #", exp_i)
  truth_df <- populateDF(WETA_filtered, covObj$siteCovs, covObj$obsCovs, unique(WETA_filtered$site), TRUE_OCC_COEFF, TRUE_DET_COEFF)  
  
  # WETA_2017$site <- -1
  # WETA_filtered only has site assignments, not the actual
  # changed occ variables; confirmed and changed 4/22
  
  # TODO: 
  #   * WETA_filtered unique locations
  #   * change clusterStats to link Checklists before running thru occ model
  
  WETA_filtered$vertex <- seq(1:nrow(WETA_filtered))
  
  og_data <- subset(WETA_filtered, select = -c(site))
  
  og_data <- og_data[order(og_data$checklist_id),]
  truth_df <- truth_df[order(truth_df$checklist_id),]
  
  og_data$species_observed_syn <- truth_df$species_observed_syn
  stopifnot(og_data$checklist_id == truth_df$checklist_id)
  if(exp_i == 1){
      
    list_of_clusterings <- baselineExp(tests, og_data, covObj, truth_df)
    

    #######
    # CREATE CONCLUSTERS
    #######
    WETA_sites <- subset(og_data, select = c(checklist_id, vertex))

    df_to_join <- list()
    tests_to_group <- c("DBSC", "rounded-4", "eBird_simple")
    for(j in 1:length(list_of_clusterings)){
      if(names(list_of_clusterings[j]) %in% tests_to_group){
        df_to_join <- append(df_to_join, list_of_clusterings[j])
      }
    }

    k <- 1
    for(df in df_to_join){
      #changing join on checklist_id to vertex because checklist_id may no longer be unique (boot exp)
      WETA_sites <- inner_join(WETA_sites, df[c("vertex", "site")], by="vertex", suffix=c("", as.character(k)))
      k <- k + 1
    }

    #TODO: 4/13 speed up?
    ptm <- proc.time()
    comb_df <- combineMethods(proj_cent, WETA_sites, og_data)
    end <- proc.time() - ptm
    disp("BALLS alg duration: ", as.character(end))

    ptm <- proc.time()
    comb_aggl_df_fast <- combineMethodsAgg(proj_cent, WETA_sites, og_data, run_mod = TRUE)
    end <- proc.time() - ptm
    disp("AGGLOM-UPDATED alg duration: ", as.character(end))
    #######

    list_of_clusterings[[paste0("balls-", paste0(tests_to_group, collapse = '-'))]] <- comb_df
    list_of_clusterings[[paste0("agglom-", paste0(tests_to_group, collapse = '-'))]] <- comb_aggl_df_fast


    #######
    # CREATE CONCLUSTERS
    #######
    WETA_sites <- subset(og_data, select = c(checklist_id, vertex))

    df_to_join <- list()
    tests_to_group <- c("clustGeo-.8-850", "DBSC", "eBird_simple")
    for(j in 1:length(list_of_clusterings)){
      if(names(list_of_clusterings[j]) %in% tests_to_group){
        df_to_join <- append(df_to_join, list_of_clusterings[j])
      }
    }

    k <- 1
    for(df in df_to_join){
      #changing join on checklist_id to vertex because checklist_id may no longer be unique (boot exp)
      WETA_sites <- inner_join(WETA_sites, df[c("vertex", "site")], by="vertex", suffix=c("", as.character(k)))
      k <- k + 1
    }

    #TODO: 4/13 speed up?
    ptm <- proc.time()
    comb_df <- combineMethods(proj_cent, WETA_sites, og_data)
    end <- proc.time() - ptm
    disp("BALLS alg duration: ", as.character(end))

    ptm <- proc.time()
    comb_aggl_df_fast <- combineMethodsAgg(proj_cent, WETA_sites, og_data, run_mod = TRUE)
    end <- proc.time() - ptm
    disp("AGGLOM-UPDATED alg duration: ", as.character(end))
    #######

    list_of_clusterings[[paste0("balls-", paste0(tests_to_group, collapse = '-'))]] <- comb_df
    list_of_clusterings[[paste0("agglom-", paste0(tests_to_group, collapse = '-'))]] <- comb_aggl_df_fast


    #######
    # CREATE CONCLUSTERS
    #######
    WETA_sites <- subset(og_data, select = c(checklist_id, vertex))

    df_to_join <- list()
    tests_to_group <- c("kmSq-1000", "DBSC", "eBird_simple")
    for(j in 1:length(list_of_clusterings)){
      if(names(list_of_clusterings[j]) %in% tests_to_group){
        df_to_join <- append(df_to_join, list_of_clusterings[j])
      }
    }

    k <- 1
    for(df in df_to_join){
      #changing join on checklist_id to vertex because checklist_id may no longer be unique (boot exp)
      WETA_sites <- inner_join(WETA_sites, df[c("vertex", "site")], by="vertex", suffix=c("", as.character(k)))
      k <- k + 1
    }

    #TODO: 4/13 speed up?
    ptm <- proc.time()
    comb_df <- combineMethods(proj_cent, WETA_sites, og_data)
    end <- proc.time() - ptm
    disp("BALLS alg duration: ", as.character(end))

    ptm <- proc.time()
    comb_aggl_df_fast <- combineMethodsAgg(proj_cent, WETA_sites, og_data, run_mod = TRUE)
    end <- proc.time() - ptm
    disp("AGGLOM-UPDATED alg duration: ", as.character(end))
    #######

    list_of_clusterings[[paste0("balls-", paste0(tests_to_group, collapse = '-'))]] <- comb_df
    list_of_clusterings[[paste0("agglom-", paste0(tests_to_group, collapse = '-'))]] <- comb_aggl_df_fast

  }

  for(j in 1:length(list_of_clusterings)){
    clustering <- list_of_clusterings[[j]]
    t_name <- names(list_of_clusterings)[[j]]
    # t_name <- strsplit(clustering_name, "-")[[1]]
    
    if(!is.na(clustering)){
      res <- clusterStats(clustering, truth_df, og_data, t_name, covObj, TRUE_OCC_COEFF, TRUE_DET_COEFF, full_df = WETA_filtered, prob_MSE = TRUE)
      # if(is.na(results[[t_name]])){
      if(exp_i == 1){
        results[[t_name]] <- new("Clustering",
                 name=t_name,
                 occ.mse=list(res$mse$occ),
                 det.mse=list(res$mse$det),
                 prob.occ.mse=list(res$prob_mse$occ),
                 prob.det.mse=list(res$prob_mse$det),
                 checklist_df=clustering,
                 svs=res$svs,
                 acc.svs=res$acc.svs,
                 simil.GT=list(res$simil.GT))
        
      } else {
        results[[t_name]]@occ.mse <- append(results[[t_name]]@occ.mse, res$mse$occ)
        results[[t_name]]@det.mse <- append(results[[t_name]]@det.mse, res$mse$det)
        results[[t_name]]@prob.occ.mse <- append(results[[t_name]]@prob.occ.mse, res$prob_mse$occ)
        results[[t_name]]@prob.det.mse <- append(results[[t_name]]@prob.det.mse, res$prob_mse$det)
        results[[t_name]]@svs <- results[[t_name]]@svs + res$svs
        results[[t_name]]@acc.svs <- results[[t_name]]@acc.svs + res$acc.svs
        results[[t_name]]@simil.GT <- append(results[[t_name]]@simil.GT, list(res$simil.GT))
      }
    }
  }
}

saveRDS(results, "baseline_results_diff_combos_fixed_coeff.RDS")


# z <- combineIntoDF(results, numExps)
# 
# write.csv(z, paste0("baseline_results-", numExps, "runs.csv"))

results <- readRDS("aggregate clustering/results/cluster-results/baseline_results_diff_combos_fixed_coeff.RDS")
# results2 <- readRDS("aggregate clustering/results/cluster-results/baseline_results_diff_combos_fixed_coeff_60more.RDS")

res_mean <- list()
for(exp_i in names(results)){
  if(!(exp_i %in% c("true det coeff", "true occ coeff"))){
    res_mean[[as.character(exp_i)]]$occ <- paste0(mean(as.numeric(results[[exp_i]]@occ.mse)), "#", std(as.numeric(results[[exp_i]]@occ.mse)))
    res_mean[[as.character(exp_i)]]$det <- paste0(mean(as.numeric(results[[exp_i]]@det.mse)), "#", std(as.numeric(results[[exp_i]]@det.mse)))
    
    res_mean[[as.character(exp_i)]]$prob_occ <- paste0(mean(as.numeric(results[[exp_i]]@prob.occ.mse)), "#", std(as.numeric(results[[exp_i]]@prob.occ.mse)))
    res_mean[[as.character(exp_i)]]$prob_det <- paste0(mean(as.numeric(results[[exp_i]]@prob.det.mse)), "#", std(as.numeric(results[[exp_i]]@prob.det.mse))) 
    
    # res_mean[[as.character(exp_i)]]$simil.df <- results[[exp_i]]@simil.GT[[1]]
  }
  # }
}
row <- 1
for(exp_i in names(results)){
  if(!(exp_i %in% c("true det coeff", "true occ coeff"))){
    if(row == 1){
      simil.df <- results[[exp_i]]@simil.GT[[1]]
    } else {
      simil.df <- rbind(simil.df, results[[exp_i]]@simil.GT[[1]])
    }
    row <- row + 1
  }
}
write.csv(simil.df, "aggregate clustering/results/cluster-results/baseline/std/simil_5-20.csv")

# CC with {X, rounded-4, eBird_simple} w/ X in {DBSC, CG-.8, kmSq-1000}
# results in the same agglom clustering and nearly the same balls clustering


melt.results <- reshape2::melt(res_mean)
#
occ_res <- melt.results[melt.results$L2 == "occ",]
det_res <- melt.results[melt.results$L2 == "det",]
#

occ_res <- occ_res %>%
  separate(value, c("mean", "std"), "#")
#
det_res <- det_res %>%
  separate(value, c("mean", "std"), "#")
#
#
colnames(occ_res) <- c("mean", "std", "d/o", "type")
colnames(det_res) <- c("mean", "std", "d/o", "type")
#
occ_res$mean <- as.numeric(occ_res$mean)
occ_res$std <- as.numeric(occ_res$std)
#
det_res$mean <- as.numeric(det_res$mean)
det_res$std <- as.numeric(det_res$std)
#
library(ggplot2)
ggplot(occ_res,
       aes(x = type,
           y = `mean`,
           color = type, group = type)) +
  # geom_line() +
  geom_bar(position = position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean - std, ymax=mean+std )) +
  ggtitle("clustering alg vs occupancy estimate 40 runs")

# filtering bc some throw off the axis
occ_res.fil <- occ_res[!occ_res$type %in% c("eBird", "kmSq-1000", "DBSC"),]
ggplot(occ_res.fil,
       aes(x = type,
           y = `mean`,
           color = type, group = type)) +
  # geom_line() +
  geom_bar(position = position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean - std, ymax=mean+std )) +
  ggtitle("clustering alg vs occupancy estimate 40 runs")



det_res.fil <- det_res[!det_res$type %in% c("eBird", "kmSq-1000", "DBSC"),]
ggplot(det_res.fil,
       aes(x = type,
           y = `mean`,
           color = type, group = type)) +
  # geom_line() +
  geom_bar(position = position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean - std, ymax=mean+std )) +
  ggtitle("clustering alg vs detection estimate 40 runs")




ggplot(det_res,
       aes(x = type,
           y = `mean`,
           color = type, group = type)) +
  # geom_line() +
  geom_bar(position = position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean - std, ymax=mean+std )) +
  ggtitle("clustering alg vs detection estimate 40 runs")




res <- sqldf("SELECT * from occ_res join det_res on occ_res.type == det_res.type")
write.csv(res, "aggregate clustering/results/cluster-results/baseline/std/fixed_5-20_coeff.csv")




