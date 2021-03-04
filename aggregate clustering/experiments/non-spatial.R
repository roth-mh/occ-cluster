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

tests.A <- list(rounded=NA, eBird=NA, 
                eBird_simple=NA,
                DBSC=NA,
                clustGeo=list(),
                GT=NA,
                noisy_gt=NA,
                agnes=c(850),
                kmeans=c(850),
                local=list()
                #ANOTHER HERE?,
)

test_names.A <- list("balls", 
                     "agglomerative", 
                     # "eBird_simple",
                     "agnes-850",
                     "kmeans-850")

tests.B <- list(rounded=NA, eBird=NA, 
                eBird_simple=NA,
                DBSC=NA,
                clustGeo=list(c(.8,850)),
                GT=NA,
                noisy_gt=NA,
                agnes=c(850),
                kmeans=c(850),
                local=list()
                #ANOTHER HERE?,
)

test_names.B <- list("balls", 
                     "agglomerative", 
                     "clustGeo-.8-850",
                     "agnes-850",
                     "kmeans-850")

versions <- c("A","B")


####################
#### PRE LOADING OBJ
####################
preLoad <- list()
# saveRDS(preLoad, file = "preLoad-3-3.rds")
preLoad <- readRDS("preLoad-3-3.rds")
p_exp1 <- readRDS("pvsExp2-22.rds")
TRUE_OCC_COEFF <- runif(6, -1.5, 1.5)
TRUE_DET_COEFF <- runif(6, -1.5, 1.5)

length(unique(p_exp1[[3]]$site))
length(unique(p_exp1[[4]]$site))

proj_cent <- p_exp1[[3]]
WETA_filtered <- p_exp1[[4]]

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

numExps <- 10
for(i in 1:numExps){
  
  WETA_2017$site <- -1
  WETA_2017$vertex <- seq(1:nrow(WETA_2017))
  
  og_data <- subset(WETA_2017, select = -c(site))
  og_data$species_observed_syn <- truth_df$species_observed_syn
     
  WETA_sites <- subset(og_data, select = c(checklist_id, vertex))
               
  res_obj.A <- runExp(tests.A, covObj, WETA_sites, og_data, truth_df, proj_cent)
  res_obj.B <- runExp(tests.A, covObj, WETA_sites, og_data, truth_df, proj_cent)
  
  #####
  # similarity btwn aggl method
  # and GT partition
  #####
  comp_to_ground_truth_df.A <- makeSIMIL_TO_GT.DF(res_obj.A, test_names.A)
  comp_to_ground_truth_df.B <- makeSIMIL_TO_GT.DF(res_obj.B, test_names.B)
  
  #####
  # similarity btwn aggl method
  # and each input method
  #####
  simil_input_method_df.A <- makeCLUSTER_COMP.DF(res_obj.A, test_names.A)
  simil_input_method_df.B <- makeCLUSTER_COMP.DF(res_obj.B, test_names.B)
  
  #####
  # mse values
  #####
  mse_df.A <- makeMSE.DF(res_obj.A, test_names.A)
  mse_df.B <- makeMSE.DF(res_obj.B, test_names.B)
  
  exp_type <- "non-spatial" # MUST BE distinct OR overlap OR noisy OR pCorr non-spatial OR <empty>
  exp_num <- i
  exp_name <- paste0("-", exp_type, "-", exp_num)
  write.csv(comp_to_ground_truth_df.A, file=paste0("aggregate clustering/results/", exp_type, "/similarity_to_GT_Exp", exp_name, "-A.csv"))
  write.csv(simil_input_method_df.A, file=paste0("aggregate clustering/results/", exp_type, "/similarity_to_agglom_Exp", exp_name, "-A.csv"))
  write.csv(mse_df.A, file=paste0("aggregate clustering/results/", exp_type, "/mse", exp_name, "-A.csv"))
  
  write.csv(comp_to_ground_truth_df.B, file=paste0("aggregate clustering/results/", exp_type, "/similarity_to_GT_Exp", exp_name, "-B.csv"))
  write.csv(simil_input_method_df.B, file=paste0("aggregate clustering/results/", exp_type, "/similarity_to_agglom_Exp", exp_name, "-B.csv"))
  write.csv(mse_df.B, file=paste0("aggregate clustering/results/", exp_type, "/mse", exp_name, "-B.csv"))

  if(gen.new.df){
    TRUE_OCC_COEFF <- runif(6, -1.5, 1.5)
    TRUE_DET_COEFF <- runif(6, -1.5, 1.5)
    truth_df <- populateDF(WETA_filtered, covObj$siteCovs, covObj$obsCovs, unique(WETA_filtered$site), TRUE_OCC_COEFF, TRUE_DET_COEFF)      
  }
  
}





