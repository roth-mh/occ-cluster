# testing environment to see if there is 
# difference between two or more types of clustering
# given a random set of DET/OCC coefficient pairs

setwd("/Users/MarkRoth/Documents/Oregon State/Research/eBird/occ and grouping checklists/occ-cluster/")
library(hash)
source("helper/helpers.R")
source("helper/ClustGeo/clustGeoHelper.R")
source("helper/ClustGeo/clustGeoSuite.R")
source("helper/gen_sites.R")
source("run experiments/prelimStudyHelper.R")


MIN_OBS <- 1
MAX_OBS <- 100000
####################
w_o <- load.WETA()
WETA_2017 <- w_o[[1]]
covObj <- w_o[[2]]

########################
# LOAD MANUAL TRUTH DF #
########################
WETA_filtered <- load.WETA_filtered()

createFinalSummary <- function(lst_of_tests){
  final_summary <- hash()
  final_summary[['det_coeff']] <- list()
  final_summary[['occ_coeff']] <- list()
  final_summary[['DETprob']] <- list()
  final_summary[['OCCprob']] <- list()
  final_summary[['AVG.zds']] <- list()
  final_summary[['AVG.zdOccs']] <- list()
  final_summary[['AVG.pDet']] <- list()
  for(tst in lst_of_tests){
    final_summary[[as.character(tst)]] <- list()
  }
  return(final_summary)
}

##################
# function to make 
# adding to summary 
# object easier
##################
add_to_summary <- function(res_obj, final_s, TRUE_OCC_COEFF, TRUE_DET_COEFF, list_of_tests){
  final_s[['det_coeff']] <- append(final_s[['det_coeff']], list(TRUE_DET_COEFF))
  final_s[['occ_coeff']] <- append(final_s[['occ_coeff']], list(TRUE_OCC_COEFF))
  final_s[['OCCprob']] <- append(final_s[['OCCprob']], res_obj$occupied_prob)
  final_s[['DETprob']] <- append(final_s[['DETprob']], res_obj$det_prob)
  final_s[['AVG.pDet']] <- append(final_s[['AVG.pDet']], unlist(res_obj$p_det))
  final_summary[['AVG.zds']] <- append(final_summary[['AVG.zds']], unlist(res_obj$zds))
  final_summary[['AVG.zdOccs']] <- append(final_summary[['AVG.zdOccs']], res_obj$zdOccS)
  
  for(test in list_of_tests){
    final_s[[as.character(test)]] <- append(final_s[[as.character(test)]], unlist(res[[as.character(test)]]))
  }
  
  return(final_s)
}

lst_of_tests <- c('base', 'eBird', 'eBird_simple')
final_summary <- createFinalSummary(lst_of_tests)


numIters <- 1
numExps <- 1
for(j in 1:numExps){
  TRUE_OCC_COEFF <- rnorm(6)
  TRUE_DET_COEFF <- rnorm(6)
  disp("experiment #:", as.character(j))
  
  # res <- list(det_prob=list(), occupied_prob=list(), base_occ_mse=list(), clustGeo9_occ_mse=list(), clustGeo6_occ_mse=list(), zds=list(), zdOccS=list(), p_det=list())
  res <- hash()
  res[['det_prob']] <- list()
  res[['occupied_prob']] <- list()
  res[['zds']] <- list()
  res[['zdOccS']] <- list()
  res[['p_det']] <- list()
  for(test in lst_of_tests){
    res[[as.character(test)]] <- list()
  }
    
  for(i in 1:numIters){
    disp("iteration #:", as.character(i))
  
    truth_df <- populateDF(WETA_filtered, covObj$siteCovs, covObj$obsCovs, unique(WETA_filtered$site), TRUE_OCC_COEFF, TRUE_DET_COEFF)
    n_dets <- sum(truth_df$species_observed_syn)
    n_occ_sites <- sum(truth_df$occupied)
  
    p_det <- n_dets/nrow(truth_df)
  
    # disp("% detections: ", as.character(p_det))
    ########################
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
  
    # results_obj <- runIDExp(WETA_2017, truth_df, covObj, geoClustAlpha = .8, occ_coef = TRUE_OCC_COEFF, det_coeff = TRUE_DET_COEFF)
    results_obj <- runeBirdExp(WETA_2017, truth_df, covObj, occ_coef = TRUE_OCC_COEFF, det_coeff = TRUE_DET_COEFF, lst_of_tests)
  
    # THESE TWO ARE UNCHANGING 
    res[['det_prob']] <- mean(truth_df$det_prob)
    res[['occupied_prob']] <- mean(uniq$occupied_prob)
    res[['zds']] <- append(res[['zds']], zds)
    res[['zdOccS']] <- append(res[['zdOccS']], zdOccS)
    res[['p_det']] <- append(res[['p_det']], p_det)
  
    for(test in lst_of_tests){
      res[[as.character(test)]] <- append(res[[as.character(test)]], results_obj[[as.character(test)]])
    }
  }

  final_summary <- add_to_summary(res, final_summary, TRUE_OCC_COEFF, TRUE_DET_COEFF, lst_of_tests)
}


###############################
# plotting MSE vs various stats
###############################
plot(final_summary$DETprob, final_summary$avg1.MSE)
plot(final_summary$OCCprob, final_summary$avg1.MSE)
plot(final_summary$avg.pDet, final_summary$avg1.MSE)
plot(final_summary$avg.zds, final_summary$avg1.MSE)
plot(final_summary$avg.zdOccS, final_summary$avg1.MSE)


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


save(final_summary, file="final_summary.RData")
