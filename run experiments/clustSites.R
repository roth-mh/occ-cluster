source("/nfs/stak/users/rothmark/Documents/siteDef/helpers.R")
source("/nfs/stak/users/rothmark/Documents/siteDef/clustGeoHelper.R")
source("/nfs/stak/users/rothmark/Documents/siteDef/mstHelper.R")
source("/nfs/stak/users/rothmark/Documents/siteDef/SKATERHelper.R")
source("/nfs/stak/users/rothmark/Documents/siteDef/kmsq.R")
source("/nfs/stak/users/rothmark/Documents/siteDef/prelimStudyHelper.R")
source("/nfs/stak/users/rothmark/Documents/siteDef/gen_sites.R")
source("/nfs/stak/users/rothmark/Documents/siteDef/dbscan.R")
source("/nfs/stak/users/rothmark/Documents/siteDef/dbscanSuite.R")
source("/nfs/stak/users/rothmark/Documents/siteDef/clustGeoSuite.R")
source("/nfs/stak/users/rothmark/Documents/siteDef/SKATERSuite.R")
source("/nfs/stak/users/rothmark/Documents/siteDef/DBSCHelper.R")

####################
# GLOBAL VARIABLES #
####################
MIN_OBS <- 1
MAX_OBS <- 100000

TRUE_OCC_COEFF <- c(-.5, 0.9772027357609328, -0.028406966644699994, 0.946381131416163, 0.7298163568647422, 0.20654714056879497)
TRUE_DET_COEFF <- c(-1, 0.3935988905626474, 0.6282290623116853, -0.21418691151584734, 1.0118138402429864, 0.8490091884815945)

NUM_SITES <- 0
# ratio <- 4
####################

f_in_WETA <- "/nfs/stak/users/rothmark/Documents/siteDef/syn_species_1_2017.csv"
WETA_2017_all <- read.delim(f_in_WETA, header=TRUE, sep = ",")
# WETA_2017 <- groomDataframe(WETA_2017, covObj$det_cov, covObj$occ_cov, syn_spec = T)
f_in_syn_spec_form <- "/nfs/stak/users/rothmark/Documents/siteDef/syn_species_1_formula.txt"
covObj <- loadCovariates(f_in_syn_spec_form)
WETA_2017_region <- subset(WETA_2017_all, WETA_2017_all$latitude <= 44.5)
WETA_2017_region <- subset(WETA_2017_region, WETA_2017_region$longitude <= -123)
WETA_2017 <- WETA_2017_region


f_name <- "/nfs/stak/users/rothmark/Documents/siteDef/clusteredSites_2020-12-26_.csv"
clusted_sites <- read.delim(f_name, sep=",")

WETA_filtered <- WETA_2017
WETA_filtered$site <- -1
# link sites with their respective checklists
for(row in 1:nrow(WETA_filtered)){
  # disp(row)
  long <- WETA_filtered[row, ]$longitude
  lat <- WETA_filtered[row, ]$latitude
  site_obj <- clusted_sites[clusted_sites$longitude == long & clusted_sites$latitude == lat,]
  WETA_filtered[row,]$site <- as.character(site_obj$site)
}

li_1000 <- runMultExp(WETA_filtered, covObj, 10, WETA_2017, 900)

save(li_1000, file = "/nfs/stak/users/rothmark/Documents/siteDef/PENpoint01clustSite1000.Rdata")

