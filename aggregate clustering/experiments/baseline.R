####
# file to run baseline preliminary experiments
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

MDD <- 250    # maximum detection distance in meters
MIN_OBS <- 1
MAX_OBS <- 100000

obj <- load.WETA()
WETA_2017 <- obj[[1]]
covObj <- obj[[2]]

preLoad <- list()
# saveRDS(preLoad, file = "preLoad-3-3.rds")

# this object stores 2 separate dfs:
#   1. proj_cent:       checklists grouped together based on MDD (imprecise because
#                         the groups were calculated with chaining, so sites have no
#                         maximum distance)
#   2. ground_truth.df:   checklists with ground truth sites in the `site`` column
# ususally calculating (1) takes a long time which is why I saved it as an .rds, but
# it will have to be rerun with prepProj.centers() if MDD != 250
# preLoad <- readRDS("preLoad-3-3.rds")
# p_exp1 <- readRDS("pvsExp2-22.rds")

# something is wrong with Baseline experiment,,,, 
# but only sometimes

adv_dataset <- DBSC.adv.df


if(length(preLoad) == 0){
  lst <- prepProj.centers(adv_dataset, covObj, MDD, MIN_OBS, MAX_OBS)
  proj_cent <- lst[[1]]
  ground_truth.df <- lst[[2]]
} else {
  proj_cent <- preLoad[[1]]
  ground_truth.df <- preLoad[[2]]
}

covObj$siteCovs <- c("TCA", "TCB", "TCW", "TCG")

# ground_truth.df <- rounded_norm_det.df
ground_truth.df$occupied_prob <- -1
ground_truth.df$occupied <- -1
ground_truth.df$species_observed_syn <- -1
ground_truth.df$det_prob <- -1
ground_truth.df$formatted_date <- as.Date("01-01-2011")
# ground_truth.df


CG1 <- "clustGeo-.8-400"
# only 1/alg for now
test_names <- list(
  "na1",
  "na2",
  "base",
  # "clustGeo-.8-850",
  CG1,
  "eBird_simple",
  # "eBird",
  # "eBird_upper",
  # "eBird_lower"
  "rounded-4",
  "kmSq-1000",
  "DBSC"
)

# rounded 4 affects about a ~50m^2 region

tests <- genTests(test_names)
# `randomized`: means each experiment the coefficients were randomized
# `fixed`:      means each experiment had a single pairing of occ/det coefficients
exp_name <- "randomized-50-dbsc"
results <- list()
TRUE_OCC_COEFF <- runif(6, -1.5, 1.5)
TRUE_DET_COEFF <- runif(6, -1.5, 1.5)
numExps <- 50
for(exp_i in 1:numExps){
  disp("iteration #", exp_i)
  
  if(grepl("randomized", exp_name, fixed = TRUE)){
    TRUE_OCC_COEFF <- runif(6, -1.5, 1.5)
    TRUE_DET_COEFF <- runif(6, -1.5, 1.5)
  }
  
  truth_df <- populateDF(ground_truth.df, covObj$siteCovs, covObj$obsCovs, unique(ground_truth.df$site), TRUE_OCC_COEFF, TRUE_DET_COEFF)  
  print("finished with construction of ground truth df")
  # TODO: 
  #   * ground_truth.df unique locations
  #   * change clusterStats to link Checklists before running thru occ model
  
  ground_truth.df$vertex <- seq(1:nrow(ground_truth.df))
  
  og_data <- subset(ground_truth.df, select = -c(site))
  
  og_data <- og_data[order(og_data$latitude),]
  truth_df <- truth_df[order(truth_df$latitude),]
  
  og_data$species_observed_syn <- truth_df$species_observed_syn
  stopifnot(og_data$checklist_id == truth_df$checklist_id)
  if(exp_i == 1){
      
    list_of_clusterings <- baselineExp(tests, og_data, covObj, truth_df)
    

    #######
    # CREATE CONCLUSTERS
    #######
    WETA_sites <- subset(og_data, select = c(checklist_id, vertex))

    df_to_join <- list()
    for(i in 1:3){
      if(i == 1){
        tests_to_group <- c("DBSC", "rounded-4", "eBird_simple")  
      } else if (i == 2){
        tests_to_group <- c("DBSC", "rounded-4", "clustGeo-.8-850")  
      } else if (i == 3){
        tests_to_group <- c("DBSC", "rounded-4", "clustGeo-.8-850", "kmSq-1000")  
      }
      
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
    
  }

  for(j in 1:length(list_of_clusterings)){
    clustering <- list_of_clusterings[[j]]
    t_name <- names(list_of_clusterings)[[j]]
    
    if(!is.na(clustering)){
      res <- clusterStats(clustering, truth_df, og_data, t_name, covObj, TRUE_OCC_COEFF, TRUE_DET_COEFF, full_df = ground_truth.df, prob_MSE = TRUE)
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

# saveRDS(results, "baseline_TassCap_vars.RDS")

# results <- readRDS("aggregate clustering/results/cluster-results/baseline_results_diff_combos_fixed_coeff.RDS")

res_mean <- list()
for(exp_i in names(results)){
  if(!(exp_i %in% c("true det coeff", "true occ coeff"))){
    res_mean[[as.character(exp_i)]]$occ <- paste0(mean(as.numeric(results[[exp_i]]@occ.mse)), "#", std(as.numeric(results[[exp_i]]@occ.mse)))
    res_mean[[as.character(exp_i)]]$det <- paste0(mean(as.numeric(results[[exp_i]]@det.mse)), "#", std(as.numeric(results[[exp_i]]@det.mse)))

    res_mean[[as.character(exp_i)]]$prob_occ <- paste0(mean(as.numeric(results[[exp_i]]@prob.occ.mse)), "#", std(as.numeric(results[[exp_i]]@prob.occ.mse)))
    res_mean[[as.character(exp_i)]]$prob_det <- paste0(mean(as.numeric(results[[exp_i]]@prob.det.mse)), "#", std(as.numeric(results[[exp_i]]@prob.det.mse)))
  }
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
write.csv(simil.df, paste0("run experiments/adversarial/results/similarity-", exp_name, ".csv"))


melt.results <- reshape2::melt(res_mean)
#
for(i in 1:2){
  if(i == 1){
    eval_type = "basic"
    occ_res <- melt.results[melt.results$L2 == "occ",]
    det_res <- melt.results[melt.results$L2 == "det",]
  } else {
    eval_type = "probabilities"
    occ_res <- melt.results[melt.results$L2 == "prob_occ",]
    det_res <- melt.results[melt.results$L2 == "prob_det",]
  }
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
  # #
  # library(ggplot2)
  # ggplot(occ_res,
  #        aes(x = type,
  #            y = `mean`,
  #            color = type, group = type)) +
  #   # geom_line() +
  #   geom_bar(position = position_dodge(), stat="identity") +
  #   geom_errorbar(aes(ymin=mean - std, ymax=mean+std )) +
  #   ggtitle("clustering alg vs occupancy estimate 40 runs")
  # 
  # # filtering bc some throw off the axis
  # occ_res.fil <- occ_res[!occ_res$type %in% c("eBird", "kmSq-1000", "DBSC"),]
  # ggplot(occ_res.fil,
  #        aes(x = type,
  #            y = `mean`,
  #            color = type, group = type)) +
  #   # geom_line() +
  #   geom_bar(position = position_dodge(), stat="identity") +
  #   geom_errorbar(aes(ymin=mean - std, ymax=mean+std )) +
  #   ggtitle("clustering alg vs occupancy estimate 40 runs")
  # 
  # 
  # 
  # det_res.fil <- det_res[!det_res$type %in% c("eBird", "kmSq-1000", "DBSC"),]
  # ggplot(det_res.fil,
  #        aes(x = type,
  #            y = `mean`,
  #            color = type, group = type)) +
  #   # geom_line() +
  #   geom_bar(position = position_dodge(), stat="identity") +
  #   geom_errorbar(aes(ymin=mean - std, ymax=mean+std )) +
  #   ggtitle("clustering alg vs detection estimate 40 runs")
  # 
  # 
  # 
  # 
  # ggplot(det_res,
  #        aes(x = type,
  #            y = `mean`,
  #            color = type, group = type)) +
  #   # geom_line() +
  #   geom_bar(position = position_dodge(), stat="identity") +
  #   geom_errorbar(aes(ymin=mean - std, ymax=mean+std )) +
  #   ggtitle("clustering alg vs detection estimate 40 runs")
  # 
  # 
  # 
  # 
  res <- sqldf("SELECT * from occ_res join det_res on occ_res.type == det_res.type")
  write.csv(res, paste0("run experiments/adversarial/results/", exp_name, "-", eval_type, ".csv"))
# 
# 
# 
# 
}
