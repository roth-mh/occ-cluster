setwd("/Users/MarkRoth/Documents/Oregon State/Research/")
source("eBird/site def/helper/helpers.R")
source("eBird/site def/helper/clustGeoHelper.R")
source("eBird/site def/helper/mstHelper.R")
source("eBird/site def/helper/SKATERHelper.R")
source("eBird/site def/helper/kmsq.R")
source("eBird/site def/site defns/prelimStudyHelper.R")
source("eBird/site def/helper/gen_sites.R")
source("eBird/site def/helper/clustGeoSuite.R")
source("eBird/site def/helper/SKATERSuite.R")
source("eBird/site def/helper/dbscan.R")
source("eBird/site def/helper/dbscanSuite.R")
source("eBird/site def/helper/DBSCHelper.R")
# getwd()

####################
# GLOBAL VARIABLES #
####################
MIN_OBS <- 1
MAX_OBS <- 100000

# TRUE_OCC_COEFF <- c(-.5, 0.4650489593471517, 0.25523159396930695, 0.8861590762863591, -0.20144247919650793, 0.3855035300948283)
# TRUE_DET_COEFF <- c(-1, -0.3605617900907053, 0.8502497538703003, -0.012596564563653434, 0.10894040227556012, 0.4170190133546442)

NUM_SITES <- 0
# ratio <- 4
####################

TRUE_OCC_COEFF <- c(-.5, 0.9772027357609328, -0.028406966644699994, 0.946381131416163, 0.7298163568647422, 0.20654714056879497)
TRUE_DET_COEFF <- c(-.3, 0.3935988905626474, 0.6282290623116853, -0.21418691151584734, 1.0118138402429864, 0.8490091884815945)

f_in_WETA <- "eBird/Class Imbalance/generate syn spec/data/linear/syn_species_1_2017.csv"
WETA_2017 <- read.delim(f_in_WETA, header=TRUE, sep = ",")
# WETA_2017 <- groomDataframe(WETA_2017, covObj$det_cov, covObj$occ_cov, syn_spec = T)
f_in_syn_spec_form <- "eBird/Class Imbalance/generate syn spec/data/linear/syn_species_1_formula.txt"
covObj <- loadCovariates(f_in_syn_spec_form)
WETA_2017_region <- subset(WETA_2017, WETA_2017$latitude <= 44.5)
WETA_2017_region <- subset(WETA_2017_region, WETA_2017_region$longitude <= -123)
WETA_2017 <- WETA_2017_region

# WETA_2017 <- WETA_2017_region
plot(x = WETA_2017$longitude, y = WETA_2017$latitude, main = "2017 eBird Checklists of Oregon")
points(x=WETA_2017_region$longitude, y=WETA_2017_region$latitude, col="blue")
legend("bottomright", legend=c("all checklists", "checklists used in calculations"), col=c("black", "blue"), pch=1)

kmSq <- kmsq.MSE(WETA_df = WETA_2017, rad_m = 1000, filter = FALSE)
clust_objs <- runIDExp(kmSq, covObj, WETA_2017, 900)

# DBSC.df <- runDBSC(WETA_2017, covObj = covObj)
# clust_objs <- runIDExp(DBSC.df, covObj, WETA_2017, 900)
# li_500 <- runMultExp(kmSq, covObj, 1, WETA_2017, 900)



# load(file = "eBird/site def/data from SLURM/clustSite1000.Rdata")
# 
# i <- 0
# for(item in li_1000){
#   i <- i + 1
#   row <- as.data.frame(item)
#   if(i == 1){
#     out_df = row
#   } else {
#     out_df = rbind(out_df, row)
#   }
# }
# write.csv(clust_objs$c_100000$clust900$df, file = "eBird/site def/OVEREST.csv")
# write.csv(clust_objs$pvs$clust900$df, file = "eBird/site def/normEST.csv")
