########
# Aggregate Clustering based on this paper;
# https://dl.acm.org/doi/10.1145/1217299.1217303
########
# chunking based on MDD
########
library(sf)
library(sp)
library(raster)
# library(dbscan)
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

TRUE_OCC_COEFF <- c(-.5, .85, 1, .2, -.5, -1)
TRUE_DET_COEFF <- c(1, 1, -.5, 1, -1, .5)

tests <- list(rounded=NA, eBird=NA, eBird_simple=T,DBSC=NA,
              clustGeo=list(),
              GT=NA,
              noisy_gt=NA,
              agnes=list(),
              kmeans=list(),
              local=list()
              )

test_names <- list("balls", 
                   "agglomerative", 
                   "eBird_simple")

# saveRDS(p_exp, file="pvsExp3-1-ALL.rds")
# saveRDS(p_exp2, file="pvsExp2-22.rds")
p_exp <- readRDS(file = "pvsExp3-1-ALL.rds")
p_exp2 <- readRDS(file = "pvsExp2-22.rds")

# p_exp[[3]] <- sf::st_drop_geometry(p_exp[[3]])
# p_exp2[[3]] <- sf::st_drop_geometry(p_exp2[[3]])

View(proj_center_no_geo)

res_obj <- runExp(tests, covObj = covObj, pvsExp = p_exp2)
res_obj <- runExp(tests, covObj = covObj, pvsExp = p_exp_TEST)
p_exp_TEST <- p_exp
p_exp_TEST[[3]] <- p_exp_TEST[[3]][1:934,]
nrow(p_exp_TEST[[3]])
nrow(p_exp2[[3]])

######
######
# Stopping place: HPC run fails on v_by_s_df; line 622 on aggClustHelper
# its tough to know if the files are synced up bc i cannot test locally
# bc of a stack error that I might have to figure out.
######
######



# p_exp <- prepWETA.sites(WETA_2017, covObj, MDD = MDD, MIN_OBS = MIN_OBS, MAX_OBS = MAX_OBS)  

saveRDS(p_exp, file = "pvsExp3-1-ALL.rds")

i <- 1
for(df in res_obj$mse.list){
  stats_list <- calcStats(df$checklists$site, res_obj$base.mse$checklists$site)
  m_df <- as.data.frame(stats_list, row.names = test_names[i])
  i <- i + 1

  if(i == 1){
    res_df <- m_df
  } else {
    res_df <- rbind(res_df, m_df)
  }
}

res_df
# rm(res_obj)
res_obj$round$mse
res_obj$cg$mse
res_obj$balls$mse
res_obj$aggl$mse
res_obj$base$mse

WETA_sites <- res_obj$pvsExp[[1]]
og_data <- res_obj$pvsExp[[2]]


optics(og_data)





