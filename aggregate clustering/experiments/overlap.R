####
# file to test incorp. of non-spatial clustering algs
# (also test the overlap vs. new information balance)
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
                     "clustGeo-.5-850",
                     "clustGeo-.3-850"
                     )

test_names.B <- list("balls",
                     "agglom-updated",
                     "clustGeo-1-850",
                     "clustGeo-.5-850",
                     "clustGeo-0-850"
                     )

test_names.C <- list("balls",
                     "agglom-updated",
                     "clustGeo-.9-850",
                     "clustGeo-.8-850",
                     "clustGeo-.7-850"
                     )

tests.A <- genTests(test_names.A)
tests.B <- genTests(test_names.B)
tests.C <- genTests(test_names.C)


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

truth_df <- populateDF(WETA_filtered, covObj$siteCovs, covObj$obsCovs, unique(WETA_filtered$site), TRUE_OCC_COEFF, TRUE_DET_COEFF)  
gen.new.df <- TRUE
local_test <- FALSE

list_of_res.A <- vector("list", length(test_names.A))
list_of_res.B <- vector("list", length(test_names.B))
list_of_res.C <- vector("list", length(test_names.B))

aggl_df.A <- NA
balls_df.A <- NA

aggl_df.B <- NA
balls_df.B <- NA

aggl_df.C <- NA
balls_df.C <- NA

numExps <- 1
for(i in 1:numExps){
  
  similarity_between_inputs_algs_df.A <- NA
  similarity_between_inputs_algs_df.B <- NA
  similarity_between_inputs_algs_df.C <- NA
  
  disp("running experiment #:", as.character(i))
  
  WETA_2017$site <- -1
  WETA_2017$vertex <- seq(1:nrow(WETA_2017))
  
  og_data <- subset(WETA_2017, select = -c(site))
  
  og_data <- og_data[order(og_data$checklist_id),]
  truth_df <- truth_df[order(truth_df$checklist_id),]
  
  og_data$species_observed_syn <- truth_df$species_observed_syn
  
  WETA_sites <- subset(og_data, select = c(checklist_id, vertex))
  
  if(local_test){
    samp <- sample(seq(1:nrow(truth_df)), 200)
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
  res_obj.C <- runExp(tests.C, covObj, W_sites, og_d, t_df, p_cent, comb_df = balls_df.C, comb_aggl_df_fast = aggl_df.C)
  
  aggl_df.A <- res_obj.A$aggl_df
  balls_df.A <- res_obj.A$balls_df
  
  aggl_df.B <- res_obj.B$aggl_df
  balls_df.B <- res_obj.B$balls_df
  
  aggl_df.C <- res_obj.C$aggl_df
  balls_df.C <- res_obj.C$balls_df
  #####
  # average similarity btwn every
  # input alg
  #####
  input_simil.A <- makeINPUT_SIMUL.DF(res_obj.A)
  input_simil.B <- makeINPUT_SIMUL.DF(res_obj.B)
  input_simil.C <- makeINPUT_SIMUL.DF(res_obj.C)
  
  
  #####
  # similarity btwn aggl method
  # and GT partition
  #####
  comp_to_ground_truth_df_i.A <- makeSIMIL_TO_GT.DF(res_obj.A, test_names.A)
  comp_to_ground_truth_df_i.B <- makeSIMIL_TO_GT.DF(res_obj.B, test_names.B)
  comp_to_ground_truth_df_i.C <- makeSIMIL_TO_GT.DF(res_obj.C, test_names.C)
  
  #####
  # similarity btwn aggl method
  # and each input method
  #####
  # simil_input_method_df.A <- makeCLUSTER_COMP.DF(res_obj.A, test_names.A)
  # simil_input_method_df.B <- makeCLUSTER_COMP.DF(res_obj.B, test_names.B)
  
  #####
  # mse values
  #####
  mse_df_i.A <- makeMSE.DF(res_obj.A, test_names.A)
  mse_df_i.B <- makeMSE.DF(res_obj.B, test_names.B)
  mse_df_i.C <- makeMSE.DF(res_obj.C, test_names.C)
  
  
  for(t in 1:length(test_names.A)){
    t_name <- test_names.A[t]
    t_res <- comp_to_ground_truth_df_i.A[as.character(t_name),]
    if(is.null(list_of_res.A[[t]])){
      list_of_res.A[[t]] <- t_res
    } else {
      list_of_res.A[[t]] <- rbind(list_of_res.A[[t]], t_res)
    }
  }
  
  for(t in 1:length(test_names.B)){
    t_name <- test_names.B[t]
    t_res <- comp_to_ground_truth_df_i.B[as.character(t_name),]
    if(is.null(list_of_res.B[[t]])){
      list_of_res.B[[t]] <- t_res
    } else {
      list_of_res.B[[t]] <- rbind(list_of_res.B[[t]], t_res)
    }
  }
  
    for(t in 1:length(test_names.C)){
    t_name <- test_names.C[t]
    t_res <- comp_to_ground_truth_df_i.C[as.character(t_name),]
    if(is.null(list_of_res.C[[t]])){
      list_of_res.C[[t]] <- t_res
    } else {
      list_of_res.C[[t]] <- rbind(list_of_res.C[[t]], t_res)
    }
  }
  
  
  if(i != 1){
    temp <- cbind(comp_to_ground_truth_df.A, comp_to_ground_truth_df_i.A)
    comp_to_ground_truth_df.A <- sapply(unique(colnames(temp)), function(x) rowSums(temp[, colnames(temp) == x, drop = FALSE]))  
    
    temp <- cbind(mse_df.A, mse_df_i.A)
    mse_df.A <- sapply(unique(colnames(temp)), function(x) rowSums(temp[, colnames(temp) == x, drop = FALSE]))  
    
    temp <- cbind(comp_to_ground_truth_df.B, comp_to_ground_truth_df_i.B)
    comp_to_ground_truth_df.B <- sapply(unique(colnames(temp)), function(x) rowSums(temp[, colnames(temp) == x, drop = FALSE]))

    temp <- cbind(mse_df.B, mse_df_i.B)
    mse_df.B <- sapply(unique(colnames(temp)), function(x) rowSums(temp[, colnames(temp) == x, drop = FALSE]))
    
    temp <- cbind(comp_to_ground_truth_df.C, comp_to_ground_truth_df_i.C)
    comp_to_ground_truth_df.C <- sapply(unique(colnames(temp)), function(x) rowSums(temp[, colnames(temp) == x, drop = FALSE]))

    temp <- cbind(mse_df.C, mse_df_i.C)
    mse_df.C <- sapply(unique(colnames(temp)), function(x) rowSums(temp[, colnames(temp) == x, drop = FALSE]))
  } else {
    comp_to_ground_truth_df.A <- comp_to_ground_truth_df_i.A
    mse_df.A <- mse_df_i.A
    
    comp_to_ground_truth_df.B <- comp_to_ground_truth_df_i.B
    mse_df.B <- mse_df_i.B
    
    comp_to_ground_truth_df.C <- comp_to_ground_truth_df_i.C
    mse_df.C <- mse_df_i.C
  }

  if(gen.new.df){
    TRUE_OCC_COEFF <- runif(6, -1.5, 1.5)
    TRUE_DET_COEFF <- runif(6, -1.5, 1.5)
    truth_df <- populateDF(WETA_filtered, covObj$siteCovs, covObj$obsCovs, unique(WETA_filtered$site), TRUE_OCC_COEFF, TRUE_DET_COEFF)      
  }
  
}

comp_to_ground_truth_df.A <- comp_to_ground_truth_df.A/numExps
mse_df.A <- mse_df.A/numExps

comp_to_ground_truth_df.B <- comp_to_ground_truth_df.B/numExps
mse_df.B <- mse_df.B/numExps

comp_to_ground_truth_df.C <- comp_to_ground_truth_df.C/numExps
mse_df.C <- mse_df.C/numExps

exp_type <- "overlap" # MUST BE distinct OR overlap OR noisy OR overlap non-spatial OR <empty>
exp_name <- paste0("-", exp_type, "-totals")
write.csv(comp_to_ground_truth_df.A, 
          file=paste0("aggregate clustering/results/", 
                      exp_type, 
                      "/similarity_to_GT_Exp", 
                      exp_name, 
                      "-A.csv")
          )

write.csv(comp_to_ground_truth_df.B, 
          file=paste0("aggregate clustering/results/", 
                      exp_type, 
                      "/similarity_to_GT_Exp", 
                      exp_name, 
                      "-B.csv")
)

write.csv(comp_to_ground_truth_df.C, 
          file=paste0("aggregate clustering/results/", 
                      exp_type, 
                      "/similarity_to_GT_Exp", 
                      exp_name, 
                      "-C.csv")
)

write.csv(input_simil.A, file=paste0("aggregate clustering/results/", exp_type, "/input-simil", exp_name, "-A.csv"))
write.csv(input_simil.B, file=paste0("aggregate clustering/results/", exp_type, "/input-simil", exp_name, "-B.csv"))
write.csv(input_simil.C, file=paste0("aggregate clustering/results/", exp_type, "/input-simil", exp_name, "-C.csv"))

write.csv(mse_df.A, file=paste0("aggregate clustering/results/", exp_type, "/mse", exp_name, "-A.csv"))
write.csv(mse_df.B, file=paste0("aggregate clustering/results/", exp_type, "/mse", exp_name, "-B.csv"))
write.csv(mse_df.C, file=paste0("aggregate clustering/results/", exp_type, "/mse", exp_name, "-C.csv"))
for(t in 1:length(list_of_res.A)){
  # exp_num <- i
  exp_name <- paste0("-", test_names.A[t])
  write.csv(list_of_res.A[[t]], 
            file=paste0("aggregate clustering/results/", 
                        exp_type, 
                        "/individual/similarity_to_GT_Exp", 
                        exp_name, 
                        "-A.csv")
            )
}


for(t in 1:length(list_of_res.B)){
  # exp_num <- i
  exp_name <- paste0("-", test_names.B[t])
  write.csv(list_of_res.B[[t]], 
            file=paste0("aggregate clustering/results/", 
                        exp_type, 
                        "/individual/similarity_to_GT_Exp", 
                        exp_name, 
                        "-B.csv")
  )
}

for(t in 1:length(list_of_res.C)){
  # exp_num <- i
  exp_name <- paste0("-", test_names.C[t])
  write.csv(list_of_res.C[[t]], 
            file=paste0("aggregate clustering/results/", 
                        exp_type, 
                        "/individual/similarity_to_GT_Exp", 
                        exp_name, 
                        "-C.csv")
  )
}