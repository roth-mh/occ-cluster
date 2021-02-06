# testing environment to create a static synthetic species
setwd("Documents/Oregon State/Research/eBird/occ and grouping checklists/occ-cluster/MSETestEnv/")
source("helper/helpers.R")
source("helper/clustGeoHelper.R")
source("helper/clustGeoSuite.R")
source("helper/gen_sites.R")
source("helper/prelimStudyHelper.R")
library(ggplot2)


MIN_OBS <- 1
MAX_OBS <- 100000
####################

f_in_WETA <- "2017_UPDATED_COVS_df.csv"
WETA_2017_all <- read.delim(f_in_WETA, header=TRUE, sep = ",")
f_in_syn_spec_form <- "syn_species_1_formula.txt"
covObj <- loadCovariates(f_in_syn_spec_form)
covObj$siteCovs <- as.character(c("fall_nbr_TCA_mean_75",
                                  "fall_nbr_B4_stdDev_150",
                                  "elevation_stdDev_150",
                                  "spring_nbr_B7_stdDev_300",
                                  "aspect_mean_300"))
WETA_2017_region <- subset(WETA_2017_all, WETA_2017_all$latitude <= 44.5)
WETA_2017_region <- subset(WETA_2017_region, WETA_2017_region$longitude <= -123)
WETA_2017 <- WETA_2017_region


########################
# LOAD GROUND TRUTH DF #
########################
f_name <- "clusteredSites_2020-12-26_.csv"
clusted_sites <- read.delim(f_name, sep=",")
WETA_filtered <- WETA_2017
WETA_filtered$site <- -1
WETA_filtered$det_prob <- -1

# link checklists to their respective sites
for(row in 1:nrow(WETA_filtered)){
  # disp(row)
  long <- WETA_filtered[row, ]$longitude
  lat <- WETA_filtered[row, ]$latitude
  site_obj <- clusted_sites[clusted_sites$longitude == long & clusted_sites$latitude == lat,]
  WETA_filtered[row,]$site <- as.character(site_obj$site)
}


# TRUE_OCC_COEFF <- c(.5, .85, 1, .2, -.5, -1)
# TRUE_DET_COEFF <- c(1.5, 1, -.5, 1, -1, .5)

#################################
# initialize final_summary object
#################################
final_summary <- list(det_coeff=list(), occ_coeff=list(),
                      DETprob=list(), OCCprob=list(),
                      avg.MSE=list(), avg.zds=list(),
                      avg.zdOccS=list(), avg.pDet=list())


#######
# number of iterations PER coefficient pairing
#######
numIters <- 4

#######
# number of different coefficient pairings
#######
numExps <- 5
for(j in 1:numExps) {
  disp("experiment #:", as.character(j))

  TRUE_OCC_COEFF <- rnorm(6)
  TRUE_DET_COEFF <- rnorm(6)
  res <- list(det_prob=list(), occupied_prob=list(), base_occ_mse=list(), clustGeo9_occ_mse=list(), clustGeo6_occ_mse=list(), zds=list(), zdOccS=list(), p_det=list())
  for(i in 1:numIters){
      disp("iteration #:", as.character(i))
    
      truth_df <- populateDF(WETA_filtered, covObj$siteCovs, covObj$obsCovs, unique(WETA_filtered$site), TRUE_OCC_COEFF, TRUE_DET_COEFF)
      n_dets <- sum(truth_df$species_observed_syn)
      n_occ_sites <- sum(truth_df$occupied)
      p_det <- n_dets/nrow(truth_df) 
      ########################
      # display stats
      ########################
      # disp("% detections: ", as.character(p_det))
      # hist(truth_df$det_prob)
      # hist(uniq$occupied_prob)
      # mean(truth_df$det_prob)
      # mean(uniq$occupied_prob)
      # table(uniq$occupied)
      # table(truth_df$species_observed_syn)
      ########################
      
      uniq <- sqldf("SELECT occupied_prob, occupied, site FROM truth_df GROUP BY site")
      
      uniqDETS <- sqldf("SELECT site, sum(species_observed_syn) FROM truth_df WHERE occupied = 1 GROUP BY site")
      zds <- calc_zero_det_sites(list(truth_df))[[1]]
      zdOccS <- nrow(uniqDETS[uniqDETS$`sum(species_observed_syn)` == 0,])/nrow(uniqDETS)
    
      results_obj <- runIDExp(WETA_2017, truth_df, covObj, geoClustAlpha = .8, occ_coef = TRUE_OCC_COEFF, det_coeff = TRUE_DET_COEFF)
    
      # THESE TWO ONLY CHANGE WHEN THE COEFF CHANGE
      res$det_prob <- mean(truth_df$det_prob)
      res$occupied_prob <- mean(uniq$occupied_prob)
        
        
      res$base_occ_mse <- append(res$base_occ_mse, results_obj$base.MSE)
      res$clustGeo9_occ_mse <- append(res$clustGeo9_occ_mse, results_obj$cg9.MSE)
      res$clustGeo6_occ_mse <- append(res$clustGeo6_occ_mse, results_obj$cg6.MSE)
      res$zds <- append(res$zds, zds)
      res$zdOccS <- append(res$zdOccS, zdOccS)
      res$p_det <- append(res$p_det, p_det)
    
    }

  final_summary <- add_to_summary(res, final_summary, TRUE_OCC_COEFF, TRUE_DET_COEFF)

}

###########################
# plot MSE vs various stats
###########################
plot(final_summary$DETprob, final_summary$avg.MSE)
plot(final_summary$OCCprob, final_summary$avg.MSE)
plot(final_summary$avg.pDet, final_summary$avg.MSE)
plot(final_summary$avg.zds, final_summary$avg.MSE)
plot(final_summary$avg.zdOccS, final_summary$avg.MSE)


###############################
# turning final_summary into DF
# for use in ggplot vizs
###############################
m <- matrix(nrow = length(final_summary$DETprob), ncol = length(final_summary)-2)
df <- as.data.frame(m)
colnames(df) <- names(final_summary[3:length(final_summary)])

df$DETprob <- final_summary$DETprob
df$OCCprob <- final_summary$OCCprob
df$avg.pDet <- final_summary$avg.pDet
df$avg.zds <- final_summary$avg.zds
df$avg.zdOccS <- final_summary$avg.zdOccS
df$avg.MSE <- final_summary$avg.MSE
#################################


###########################
# plotting MSE as a func of 
# detProb and occProb 
###########################
ggplot(df, aes(x=unlist(DETprob), y=unlist(OCCprob), size = unlist(avg.MSE))) +
  geom_point(alpha=0.7)




