setwd("/Users/MarkRoth/Documents/Oregon State/Research/eBird/occ and grouping checklists/occ-cluster/")
source("helper/helpers.R")
source("helper/ClustGeo/clustGeoHelper.R")
source("helper/ClustGeo/clustGeoSuite.R")
source("helper/SKATER/mstHelper.R")
source("helper/SKATER/SKATERHelper.R")
source("helper/SKATER/SKATERSuite.R")
source("helper/gen_sites.R")
source("helper/kmsq.R")
source("run experiments/prelimStudyHelper.R")
source("helper/DBSCAN/dbscan.R")
source("helper/DBSCAN/dbscanSuite.R")
source("helper/DBSC/DBSCHelper.R")
# getwd()

####################
# GLOBAL VARIABLES #
####################
MIN_OBS <- 1
# MAX_OBS <- 10
MAX_OBS <- 100000

# TRUE_OCC_COEFF <- c(-.5, 0.4650489593471517, 0.25523159396930695, 0.8861590762863591, -0.20144247919650793, 0.3855035300948283)
# TRUE_DET_COEFF <- c(-1, -0.3605617900907053, 0.8502497538703003, -0.012596564563653434, 0.10894040227556012, 0.4170190133546442)

NUM_SITES <- 0
# ratio <- 4
####################

f_in_WETA <- "../../../ICB General/data generation/2017_UPDATED_COVS_df.csv"
WETA_2017_all <- read.delim(f_in_WETA, header=TRUE, sep = ",")
# WETA_2017 <- groomDataframe(WETA_2017, covObj$det_cov, covObj$occ_cov, syn_spec = T)
f_in_syn_spec_form <- "../../Class Imbalance/generate syn spec/data/linear/syn_species_1_formula.txt"
covObj <- loadCovariates(f_in_syn_spec_form)
covObj$siteCovs <- as.character(c("fall_nbr_TCA_mean_75", 
                                  "fall_nbr_B4_stdDev_150", 
                                  "elevation_stdDev_150", 
                                  "spring_nbr_B7_stdDev_300", 
                                  "aspect_mean_300"))
WETA_2017_region <- subset(WETA_2017_all, WETA_2017_all$latitude <= 44.5)
WETA_2017_region <- subset(WETA_2017_region, WETA_2017_region$longitude <= -123)
WETA_2017 <- WETA_2017_region

# WETA_2017_all[order(as.Date(WETA_2017_all$as_date, format="%Y-%m-%d")),][nrow(WETA_2017_all),]

# WETA_2017 <- WETA_2017_region
# plot(x = WETA_2017_all$longitude, y = WETA_2017_all$latitude, main = "2017 eBird Checklists of Oregon")
# points(x=WETA_2017$longitude, y=WETA_2017$latitude, col="blue")
# legend("bottomright", legend=c("all checklists", "checklists used in calculations"), col=c("black", "blue"), pch=1)

f_name <- "../../../clusteredSites_2020-12-26_.csv"
clusted_sites <- read.delim(f_name, sep=",")
WETA_filtered <- WETA_2017
WETA_filtered$site <- -1
WETA_filtered$det_prob <- -1

# link sites with their respective checklists
for(row in 1:nrow(WETA_filtered)){
  # disp(row)
  long <- WETA_filtered[row, ]$longitude
  lat <- WETA_filtered[row, ]$latitude
  site_obj <- clusted_sites[clusted_sites$longitude == long & clusted_sites$latitude == lat,]
  WETA_filtered[row,]$site <- as.character(site_obj$site)
}

results_obj <- -1
while(results_obj == -1){
  # TRUE_OCC_COEFF <- rnorm(6)
  # TRUE_DET_COEFF <- rnorm(6)
  # I
  # TRUE_OCC_COEFF <- c(-.5, .85, 1, .2, -.5, -1) # II
  # TRUE_DET_COEFF <- c(2, 1, -.5, 1, -1, .5)
  TRUE_OCC_COEFF <- c(-.5, .85, 1, .2, -.5, -1) # II
  TRUE_DET_COEFF <- c(1, 1, -.5, 1, -1, .5)
  # p_det <- 0
  # while(p_det < .15){
  truth_df <- populateDF(WETA_filtered, covObj$siteCovs, covObj$obsCovs, unique(WETA_filtered$site), TRUE_OCC_COEFF, TRUE_DET_COEFF)
  n_dets <- sum(truth_df$species_observed_syn)
  n_occ_sites <- sum(truth_df$occupied)

  p_det <- n_dets/nrow(truth_df)
  
  disp("% detections: ", as.character(n_dets/nrow(truth_df)))
  ########################
  hist(truth_df$det_prob)
  uniq <- sqldf("SELECT occupied_prob, occupied, site FROM truth_df GROUP BY site")
  hist(uniq$occupied_prob)
  # mean(uniq$occupied_prob)
  # table(uniq$occupied)
  # table(truth_df$species_observed_syn)
  ########################
  
  uniqDETS <- sqldf("SELECT site, sum(species_observed_syn) FROM truth_df WHERE occupied = 1 GROUP BY site")
  zds <- calc_zero_det_sites(list(truth_df))[[1]]
  zdOccS <- nrow(uniqDETS[uniqDETS$`sum(species_observed_syn)` == 0,])/nrow(uniqDETS)
  zds
  zdOccS
  # }
  
  results_obj <- runSingleExp(truth_df, covObj, WETA_2017, 900)
}




# Thoughts:
#  - zdOccS (zeroDetectionOccupiedSites) may have something to do with the wonkiness, more than zds
#  - % detection after accounting for sites
#  - changing intercept shifts our dist.
#  - fall_nbr_TCA_mean_75 is q influential
#
#  - # of detections depends on the occupied sites ... can have det
#     prob = 1 and still very few detections
#  - something i haven't tried yet is a high detection probability 
#


i1.good <- list(results_obj, zds, zdOccS, n_dets/nrow(truth_df))
i2.ok 
i3.bad 
i4.ok 
# TRUE_OCC_COEFF <- c(-.5, .5, -1, 1, -.5, -1)
# TRUE_DET_COEFF <- c(.5, 1, -.5, 1, -1, .5)




print.MSE(results_obj)
print.MSE <- function(r_obj){
  for(i in 1:(length(r_obj)-1)){
    disp(names(r_obj)[[i]])
    disp(as.character(r_obj[[i]]$mse$occ))
    disp(as.character(r_obj[[i]]$mse$det))
    disp("~~~~~~~~~~~~~~~~~~~~~")
  }
}

colMeans(subset(WETA_filtered, select=c(covObj$siteCovs)))

var(WETA_filtered$fall_nbr_TCA_mean_75)
var(WETA_filtered$fall_nbr_B4_stdDev_150)
var(WETA_filtered$elevation_stdDev_150)
var(WETA_filtered$spring_nbr_B7_stdDev_300)
var(WETA_filtered$aspect_mean_300)


load(file = "../data from SLURM/PENpoint01clustSite1000.Rdata")

out_df <- as.data.frame(matrix())
i <- 0
for(item in results_obj){
  i <- i + 1
  item$name <- names(results_obj)[[i]]
  # if(length(item$zds) == 0){
  #   item$zds <- 0
  # }
  row <- as.data.frame(item)
  if(i == 1){
    out_df = row
  } else {
    out_df = rbind(out_df, row)
  }
}
write.csv(out_df, file = "../data from SLURM/PENpoint01clustSites.csv")
# write.csv(clust_objs$c_100000$clust900$df, file = "eBird/site def/OVEREST.csv")
# write.csv(clust_objs$pvs$clust900$df, file = "eBird/site def/normEST.csv")
