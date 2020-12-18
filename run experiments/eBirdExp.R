########
# Variations on eBird Best Practices:
#   - best practice
#   - best practice - observer constraint
#   - rounding to 4 decimals
#   - rounding to 4 decimals - observer constraint
########

setwd("/Users/MarkRoth/Documents/Oregon State/Year 1/Research/")
source("eBird/site def/helper/helpers.R")


f_in_WETA <- "eBird/Class Imbalance/generate syn spec/data/linear/syn_species_0_2017.csv"
WETA_2017 <- read.delim(f_in_WETA, header=TRUE, sep = ",")

f_in_syn_spec_form <- "eBird/Class Imbalance/generate syn spec/data/linear/syn_species_0_formula.txt"
covObj <- loadCovariates(f_in_syn_spec_form)




w <- WETA_2017
w1 <- WETA_2017
sites_ebird_filter <- filter_repeat_visits(
  w,
  min_obs = as.integer(2), 
  max_obs = as.integer(10),
  annual_closure = TRUE,
  date_var = "formatted_date",
  site_vars = c("locality_id", "observer_id")
)

sites_ebird_filter_no_OBS <- filter_repeat_visits(
  w1,
  min_obs = as.integer(2), 
  max_obs = as.integer(10),
  annual_closure = TRUE,
  date_var = "formatted_date",
  site_vars = c("locality_id")
)


WETA_2017_ROUND1 <- roundLatLong(WETA_2017, 4)
WETA_2017_ROUND2 <- roundLatLong(WETA_2017, 4)
sites_ebird_filter_ROUNDED <- filter_repeat_visits(
  WETA_2017_ROUND1,
  min_obs = as.integer(2), max_obs = as.integer(10),
  annual_closure = TRUE,
  date_var = "formatted_date",
  site_vars = c("rounded_locality_id", "observer_id")
)


sites_ebird_filter_ROUNDED_no_OBS <- filter_repeat_visits(
  WETA_2017_ROUND2,
  min_obs = as.integer(2), max_obs = as.integer(10),
  annual_closure = TRUE,
  date_var = "formatted_date",
  site_vars = c("rounded_locality_id")
)

w2017_stats <- siteStatsDt(sites_ebird_filter, covObj = covObj)
w2017_stats_no_OBS <- siteStatsDt(sites_ebird_filter_no_OBS, covObj = covObj)
w2017_3_stats <- siteStatsDt(sites_ebird_filter_ROUNDED, covObj = covObj)
w2017_3_no_OBS <- siteStatsDt(sites_ebird_filter_ROUNDED_no_OBS, covObj = covObj)



#################
# compare Occ MSE
# syn_spec_0
true_occ_coeff <- c(-.5, 0.4650489593471517, 0.25523159396930695, 0.8861590762863591, -0.20144247919650793, 0.3855035300948283)
true_det_coeff <- c(-1, -0.3605617900907053, 0.8502497538703003, -0.012596564563653434, 0.10894040227556012, 0.4170190133546442)
no_OBS_MSE <- calcOccMSE(sites_df = sites_ebird_filter_no_OBS, 
                         covariate_object = covObj,
                         true_occ_coefficients = true_occ_coeff,
                         true_det_coefficients = true_det_coeff,
                         syn_spec = TRUE
)

# debug(calcOccMSE)
eBird_MSE <- calcOccMSE(sites_df = sites_ebird_filter, 
                        covariate_object = covObj,
                        true_occ_coefficients = true_occ_coeff,
                        true_det_coefficients = true_det_coeff,
                        syn_spec = TRUE
)

no_OBS_ROUNDED_MSE <- calcOccMSE(sites_df = sites_ebird_filter_ROUNDED_no_OBS,
                                covariate_object = covObj,
                                true_occ_coefficients = true_occ_coeff,
                                true_det_coefficients = true_det_coeff,
                                syn_spec = TRUE
)

ROUNDED_MSE <- calcOccMSE(sites_df = sites_ebird_filter_ROUNDED, 
                          covariate_object = covObj,
                          true_occ_coefficients = true_occ_coeff,
                          true_det_coefficients = true_det_coeff,
                          syn_spec = TRUE
)

#################



eBird_MSE$occ
no_OBS_MSE$occ
ROUNDED_MSE$occ
no_OBS_ROUNDED_MSE$occ

eBird_MSE$det
no_OBS_MSE$det
ROUNDED_MSE$det
no_OBS_ROUNDED_MSE$det

